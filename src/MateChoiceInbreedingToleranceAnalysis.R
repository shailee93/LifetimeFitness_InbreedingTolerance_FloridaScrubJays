# Analysis

library(tidyverse)
library(lme4)
library(MASS)
library(glmmTMB)
library(knitr)
library(scales)
library(survival)
library(Matrix)
library(jtools)
library(DHARMa)
library(ggpattern)
library(patchwork)
library(ggeffects)

# setwd to directory with "src" folder that has this script in it
setwd("..")

# 1. Relatedness of potential mates with distance ----------------------------------

AllPossPairsData <- read_csv("./data/AllPossPairsData.csv")

## Load and run custom function to fit GLMMs with binomial error distributions with kinship coefficient to potential mate as the dependent variable and distance to potential mate as the fixed effect. Random effects of focal individual ID and year are also included.

source("./src/CustomFunctions/PotentialMateCoeffFn.R")

PotentialMateCoeffFn(df=AllPossPairsData)

# 2. Differences in simulated means under different random mating scenarios ----------

SimData <- read_csv("./data/SimData.csv")
SimData_displim <- read_csv("./data/SimData_displim.csv")
SimData_displim_nosexbias <- read_csv("./data/SimData_displim_nosexbias.csv")
BreederData <- read_csv("./data/BreederData.csv")

  ## add datatype column to each dataset

  SimData$datatype <- "Random mating with unconstrained dispersal"
  SimData_displim$datatype <- "Random mating with sex-biased limited dispersal"
  SimData_displim_nosexbias$datatype <- "Random mating with limited dispersal"
  BreederData$datatype <- "Observed"

  ## make combined dataset of simulated data

  comb_Sims_displim_nosexbias <- bind_rows(SimData, SimData_displim, SimData_displim_nosexbias)

## pairwise t-tests to test differences in kinship coefficients to mates under the different simulation scenarios

pairwise.t.test(x=comb_Sims_displim_nosexbias$Coeff, g=comb_Sims_displim_nosexbias$datatype, p.adjust.method = "bonf")

# 3. Observed kinship coefficient to mate compared to random mating scenarios -----------

## load function to compare observed means of kinship cofficient to mates to those simulated under different random mating scenarios

source("./src/CustomFunctions/ObsSimsFn.R")

## make a list of the datasets of the three simulated random mating scenarios

sim_datasets <- list(
  "Sims" = SimData,
  "Sims_displim" = SimData_displim,
  "Sims_displim_nosexbias" = SimData_displim_nosexbias
)
variable_labels <- c("Coeff", "Dist")
p_value_types <- c("more", "less")

for (sim_name in names(sim_datasets)) {
  for (var in variable_labels) {
    for (p_type in p_value_types) {

      col_name <- paste0(var, "_", stringr::str_to_lower(p_type))
      p_value_col <- paste0(col_name, "_p_value")

      result_df <- ObsSimsFn(
        ObsDF = BreederData,
        SimDF = sim_datasets[[sim_name]],
        cVars = c(var),
        cSexes = c("Male", "Female"),
        cPairType = c("First", "Later"),
        StartYear = 1990,
        EndYear = 2021,
        EmpPValueType = p_type
      )

      result_df <- result_df %>% dplyr::mutate(psignif = ifelse(p <= 0.05, "1_yes", "2_no"))
      result_df$psignif <- factor(result_df$psignif, levels=c("1_yes", "2_no"))

      lookup_df <- result_df %>%
        dplyr::select(PairType, Sex, psignif, p) %>%
        distinct()

      temp_joined <- sim_datasets[[sim_name]] %>%
        left_join(lookup_df, by = c("PairType", "Sex")) %>%
        rename(
          !!col_name := psignif,
          !!p_value_col := p
        )

      sim_datasets[[sim_name]] <- temp_joined
    }
  }
}
## consolidate results

Coeff_results <- map(sim_datasets, ~ .x %>%
                         group_by(Sex, PairType) %>%
                         summarize(
                           Coeff_more_p_value = unique(Coeff_more_p_value),
                           Coeff_less_p_value = unique(Coeff_less_p_value),
                           .groups = 'drop'
                         ))

## Access results
Coeff_results$Sims
Coeff_results$Sims_displim
Coeff_results$Sims_displim_nosexbias

# 4. Observed frequency of first, second, and third order kin mates compared to random mating scenarios -----------------------

### make df of proportion of first, second, and third order kin

comb_Sims_displim_nosexbias_PROPS <- comb_Sims_displim_nosexbias %>%
  mutate(kin_cat = case_when(
    Coeff >= 0.25 ~ "First",
    Coeff >= 0.125 ~ "Second",
    Coeff >= 0.0625 ~ "Third",
    TRUE ~ "non_kin"
  )) %>%
  group_by(datatype, Year, Sex, PairType, Rep) %>%
  count(kin_cat) %>% # Counts occurrences within each zone
  mutate(prop = n / sum(n)) %>%
  mutate(prop = replace_na(prop, 0)) %>%
  ungroup()

Sims_datatype_list <- split(comb_Sims_displim_nosexbias_PROPS, comb_Sims_displim_nosexbias_PROPS$datatype)

Sims_wide_PROP <- map(Sims_datatype_list, ~ .x %>%
                         pivot_wider(names_from = kin_cat,
                                     values_from = prop,
                                     id_cols = c(datatype, Year, Sex, PairType, Rep),
                                     values_fill = 0))

unique(comb_Sims_displim_nosexbias$datatype)

BreederData_PROPS_wide <- BreederData %>%
  mutate(kin_cat = case_when(
    Coeff >= 0.25 ~ "First",
    Coeff >= 0.125 ~ "Second",
    Coeff >= 0.0625 ~ "Third",
    TRUE ~ "non_kin"
  )) %>%
  group_by(Sex, Year, PairType) %>%
  count(kin_cat) %>% # Counts occurrences within each zone
  mutate(prop = n / sum(n)) %>%
  mutate(prop = replace_na(prop, 0)) %>%
  ungroup() %>%
  pivot_wider(names_from = kin_cat,
              values_from = prop,
              id_cols = c(Year, Sex, PairType),
              values_fill = 0)


# Define the variables and their corresponding xlabels
kin_labels <- c("First" = "Proportion of First-order Kin Pairs",
                     "Second" = "Proportion of Second-order Kin Pairs",
                     "Third" = "Proportion of Third-order Kin Pairs",
                     "non_kin" = "Proportion of Non-kin Pairs")
p_value_types <- c("more", "less")


Sims_displim_nosexbias_PROP <- Sims_wide_PROP[[1]]
Sims_displim_PROP <- Sims_wide_PROP[[2]]
Sims_PROP <- Sims_wide_PROP[[3]]

# Create empty list to store results
Sims_PROP_RESULTS <- list()

# Loop through variables and p-value types
for (var in names(kin_labels)) {
  for (p_type in p_value_types) {

    # Create result name
    result_name <- paste0(var, "_Sims_PROP_", stringr::str_to_title(p_type))

    # Run the function
    Sims_PROP_RESULTS[[result_name]] <- ObsSimsFn(
      ObsDF = BreederData_PROPS_wide,
      SimDF = Sims_PROP,
      cVars = c(var),
      cSexes = c("Male", "Female"),
      cPairType = c("First", "Later"),
      StartYear = 1990,
      EndYear = 2021,
      EmpPValueType = p_type)

    Sims_PROP_RESULTS[[result_name]]$psignif <- ifelse(Sims_PROP_RESULTS[[result_name]]$p <= 0.05, "1_yes", "2_no")
    Sims_PROP_RESULTS[[result_name]]$psignif <- factor(Sims_PROP_RESULTS[[result_name]]$psignif, levels=c("1_yes", "2_no"))

  }
}

# Create empty list to store results
Sims_displim_PROP_RESULTS <- list()

# Loop through variables and p-value types
for (var in names(kin_labels)) {
  for (p_type in p_value_types) {

    # Create result name
    result_name <- paste0(var, "_Sims_PROP_", stringr::str_to_title(p_type))

    # Run the function
    Sims_displim_PROP_RESULTS[[result_name]] <- ObsSimsFn(
      ObsDF = BreederData_PROPS_wide,
      SimDF = Sims_displim_PROP,
      cVars = c(var),
      cSexes = c("Male", "Female"),
      cPairType = c("First", "Later"),
      StartYear = 1990,
      EndYear = 2021,
      EmpPValueType = p_type)

    Sims_displim_PROP_RESULTS[[result_name]]$psignif <- ifelse(Sims_displim_PROP_RESULTS[[result_name]]$p <= 0.05, "1_yes", "2_no")
    Sims_displim_PROP_RESULTS[[result_name]]$psignif <- factor(Sims_displim_PROP_RESULTS[[result_name]]$psignif, levels=c("1_yes", "2_no"))

  }
}

# Create empty list to store results
Sims_displim_nosexbias_PROP_RESULTS <- list()

# Loop through variables and p-value types
for (var in names(kin_labels)) {
  for (p_type in p_value_types) {

    # Create result name
    result_name <- paste0(var, "_Sims_PROP_", stringr::str_to_title(p_type))

    # Run the function
    Sims_displim_nosexbias_PROP_RESULTS[[result_name]] <- ObsSimsFn(
      ObsDF = BreederData_PROPS_wide,
      SimDF = Sims_displim_nosexbias_PROP,
      cVars = c(var),
      cSexes = c("Male", "Female"),
      cPairType = c("First", "Later"),
      StartYear = 1990,
      EndYear = 2021,
      EmpPValueType = p_type
    )

    Sims_displim_nosexbias_PROP_RESULTS[[result_name]]$psignif <- ifelse(Sims_displim_nosexbias_PROP_RESULTS[[result_name]]$p <= 0.05, "1_yes", "2_no")
    Sims_displim_nosexbias_PROP_RESULTS[[result_name]]$psignif <- factor(Sims_displim_nosexbias_PROP_RESULTS[[result_name]]$psignif, levels=c("1_yes", "2_no"))

  }
}

## summarize results for first-order kin mates (can do the same for second and third order kin and non-kin mates)

Sims_PROP_RESULTS$First_Sims_PROP_More %>%
  group_by(Sex, PairType) %>%
  summarize(p = max(p))

Sims_displim_PROP_RESULTS$First_Sims_PROP_Less %>%
  group_by(Sex, PairType) %>%
  summarize(p = max(p))

Sims_displim_nosexbias_PROP_RESULTS$First_Sims_PROP_Less %>%
  group_by(Sex, PairType) %>%
  summarize(p = max(p))



## make Fig. 1 =========================================================

## Fig 1, left ###################################################

### summarize p value significance

sim_datasets <- map(sim_datasets, ~ .x %>%
                  mutate(psignif_coeff_comb = case_when(
                      Coeff_more == "1_yes" ~ Coeff_more,
                      Coeff_less == "1_yes" ~ Coeff_less,
                      TRUE ~ "2_no"
                    )
                  )
)

### Assign back to original dataset names
Sims <- sim_datasets[["Sims"]]
Sims_displim <- sim_datasets[["Sims_displim"]]
Sims_displim_nosexbias <- sim_datasets[["Sims_displim_nosexbias"]]

### combine simulation datasets
comb_Sims_displim_nosexbias <- bind_rows(Sims, Sims_displim, Sims_displim_nosexbias)

## calculate means and variances of simulated data
comb_Sims_displim_nosexbias_means <- comb_Sims_displim_nosexbias %>%
                                      group_by(datatype, Rep, Sex, PairType, psignif_coeff_comb) %>%
                                      summarize(mean_coeff = mean(Coeff),
                                                var_coeff = var(Coeff)) %>%
                                      ungroup()

## calculate means and variances of observed data
BreederData_means <- BreederData %>%
                  group_by(datatype, Sex, PairType) %>%
                  summarize(mean_coeff = mean(Coeff),
                            var_coeff = var(Coeff)) %>%
                  ungroup()

## combine observed with simulated data  -- means and variances
comb_BreederData_Sims_displim_nosexbias_MEANS <- bind_rows(BreederData_means, comb_Sims_displim_nosexbias_means)

## make simulated means of means df
means_of_means <- filter(comb_BreederData_Sims_displim_nosexbias_MEANS, !datatype == "Observed")  %>%
                  group_by(datatype, Sex, PairType) %>%
                  summarize(mean_mean_coeff = mean(mean_coeff)) %>%
                  ungroup()

p1 <- comb_BreederData_Sims_displim_nosexbias_MEANS %>%

  mutate(PairType = recode(PairType, First = "First pairings", Later = "Later pairings", Natal = "First pairings", Breeding = "Later pairings")) %>%

  mutate(datatype = factor(datatype, levels =c("Random mating with unconstrained dispersal", "Random mating with sex-biased limited dispersal", "Random mating with limited dispersal"))) %>%

  ### background distribution of relatedness of simulated pairs
  dplyr::filter(!datatype == "Observed") %>%
  dplyr::mutate(Sex = factor(Sex, levels = c("Female", "Male")),
                PairType = factor(PairType, levels = c("First pairings", "Later pairings"))) %>%
  group_by(Sex, PairType) %>%
  ggplot(aes(x = mean_coeff)) +
  geom_histogram(aes(alpha = psignif_coeff_comb, fill=datatype, color=datatype), bins = 30) +

  ## add vertical line of observed mean
  geom_vline(data = (dplyr::filter(comb_BreederData_Sims_displim_nosexbias_MEANS, datatype == "Observed") %>% mutate(PairType = recode(PairType, First = "First pairings", Later = "Later pairings"))),  color = "black", lwd=1.5, aes(xintercept = mean_coeff)) +

  ## add arrows of simulated means
  geom_segment(data =  means_of_means %>% mutate(PairType = recode(PairType, First = "First pairings", Later = "Later pairings")),
               aes(x = mean_mean_coeff, xend = mean_mean_coeff,
                   y = 400 * 1.5, yend = 400 * 1.2,
                   color = datatype),
               arrow = arrow(length = unit(0.3, "cm")),
              size = 1) +

  ## facet by sex and pairing type
  facet_wrap(vars(PairType, Sex), nrow=4, scales = "fixed") +

  ## labels and legends
  xlab("Relatedness to mate") +
  ylab("Count") +

  # color, transperency and linetype aesthetics
  scale_color_viridis_d(option="inferno", begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(option="inferno", begin = 0.2, end = 0.8, labels = label_wrap(19)) +
  scale_alpha_manual(values=c(0.5, 0)) +
  guides(alpha = "none", color = "none") +
  labs(fill = "Simulation type")+
  theme(legend.key.spacing.y = unit(0.1, "cm"))

p1

## Fig 1, right ###################################################

### Plot for First order kin ==========================

## add p_signif variable

# Define the simulation datasets to iterate over
sim_dataset_PROP <- list("Sims_PROP" = Sims_PROP,
                     "Sims_displim_PROP" = Sims_displim_PROP,
                     "Sims_displim_nosexbias_PROP" = Sims_displim_nosexbias_PROP)

# Process each simulation dataset
for (sim_name in names(sim_dataset_PROP)) {

  # Get the current simulation dataset
  current_sim <- sim_dataset_PROP[[sim_name]]

  # Define the result names for this simulation dataset
  result_names <- c(
    paste0("First_Sims_PROP_More"),
    paste0("First_Sims_PROP_Less"),
    paste0("Second_Sims_PROP_More"),
    paste0("Second_Sims_PROP_Less"),
    paste0("Third_Sims_PROP_More"),
    paste0("Third_Sims_PROP_Less"),
    paste0("non_kin_Sims_PROP_More"),
    paste0("non_kin_Sims_PROP_Less")

  )

  # Define corresponding column names
  col_names <- c(
    "psignif_K1_more",
    "psignif_K1_less",
    "psignif_K2_more",
    "psignif_K2_less",
    "psignif_K3_more",
    "psignif_K3_less",
    "psignif_NK_more",
    "psignif_NK_less"
  )

  # define results list
  results_list <- paste0(sim_name, "_RESULTS")
  results <- get(results_list)

  # Add psignif columns for each result
  for (i in seq_along(result_names)) {
    result_name <- result_names[i]
    col_name <- col_names[i]

    # Check if the result exists in the results list
    if (result_name %in% names(results)) {
      current_sim[[col_name]] <- results[[result_name]]$psignif[
        match(interaction(current_sim$Sex, current_sim$PairType),
              interaction(results[[result_name]]$Sex, results[[result_name]]$PairType))
      ]
    }
  }

  # Assign back to the appropriate variable name
  if (sim_name == "Sims_PROP") {
    Sims_PROP <- current_sim
  } else if (sim_name == "Sims_displim_PROP") {
    Sims_displim_PROP <- current_sim
  } else if (sim_name == "Sims_displim_nosexbias_PROP") {
    Sims_displim_nosexbias_PROP <- current_sim
  }
}

colnames(Sims_PROP)
unique(Sims_PROP$datatype)
colnames(Sims_displim_PROP)
unique(Sims_displim_PROP$datatype)
colnames(Sims_displim_nosexbias_PROP)
unique(Sims_displim_nosexbias_PROP$datatype)

# Define the simulation datasets to iterate over
sim_dataset_PROP <- list("Sims_PROP" = Sims_PROP,
                         "Sims_displim_PROP" = Sims_displim_PROP,
                         "Sims_displim_nosexbias_PROP" = Sims_displim_nosexbias_PROP)


## assign p-values
## if p_value is signif, keep that, otherwise no

# For loop to apply the same mutation to all dataframes
for (sim_name in names(sim_dataset_PROP)) {

  # Get the current simulation dataset
  current_sim <- sim_dataset_PROP[[sim_name]]

  current_sim <- current_sim %>%
    mutate(
      psignif_K1_comb = case_when(
        psignif_K1_more == "1_yes" ~ psignif_K1_more,
        psignif_K1_less == "1_yes" ~ psignif_K1_less,
        TRUE ~ "2_no"
      ),
      psignif_K2_comb = case_when(
        psignif_K2_more == "1_yes" ~ psignif_K2_more,
        psignif_K2_less == "1_yes" ~ psignif_K2_less,
        TRUE ~ "2_no"
      ),
      psignif_K3_comb = case_when(
        psignif_K3_more == "1_yes" ~ psignif_K3_more,
        psignif_K3_less == "1_yes" ~ psignif_K3_less,
        TRUE ~ "2_no"
      ),
      psignif_NK_comb = case_when(
        psignif_NK_more == "1_yes" ~ psignif_NK_more,
        psignif_NK_less == "1_yes" ~ psignif_NK_less,
        TRUE ~ "2_no"
      )
    )

  # Assign back to the appropriate variable name
  if (sim_name == "Sims_PROP") {
    Sims_PROP <- current_sim
  } else if (sim_name == "Sims_displim_PROP") {
    Sims_displim_PROP <- current_sim
  } else if (sim_name == "Sims_displim_nosexbias_PROP") {
    Sims_displim_nosexbias_PROP <- current_sim
  }
}


## combine simulation datasets
comb_Sims_displim_nosexbias_PROP <- bind_rows(Sims_PROP, Sims_displim_PROP, Sims_displim_nosexbias_PROP)
nrow(comb_Sims_displim_nosexbias) # 4100767
unique(comb_Sims_displim_nosexbias_PROP$datatype)

## calculate means and variances
comb_Sims_displim_nosexbias_means_PROP <- comb_Sims_displim_nosexbias_PROP %>%
  group_by(datatype, Rep, Sex, PairType, psignif_K1_comb, psignif_K2_comb, psignif_K3_comb, psignif_NK_comb) %>%
  summarize(mean_K1 = mean(First),
            mean_K2 = mean(Second),
            mean_K3 = mean(Third),
            mean_NK = mean(non_kin)
  ) %>%
  ungroup()

head(comb_Sims_displim_nosexbias_means_PROP)

unique(comb_Sims_displim_nosexbias_means_PROP$datatype)

## calculate means and variances of observed data
BreederData_means_PROP <- BreederData_PROPS_wide %>%
  mutate(datatype = "Observed") %>%
  group_by(datatype, Sex, PairType) %>%
  summarize(mean_K1 = mean(First),
            mean_K2 = mean(Second),
            mean_K3 = mean(Third),
            mean_NK = mean(non_kin)) %>%
  ungroup()

## combine observed with simulated data  -- means and variances
comb_BreederData_Sims_displim_nosexbias_PROPS <- bind_rows(BreederData_means_PROP, comb_Sims_displim_nosexbias_means_PROP)

## make simulated means of means df

means_of_PROPS <- filter(comb_BreederData_Sims_displim_nosexbias_PROPS, !datatype == "Observed")  %>%
  group_by(datatype, Sex, PairType) %>%
  summarize(mean_prop_K1 = mean(mean_K1),
            mean_prop_K2 = mean(mean_K2),
            mean_prop_K3 = mean(mean_K3),
            mean_prop_NK = mean(mean_NK),
            ) %>%
  ungroup()

# Define the relatedness types to loop through
relatedness_types <- c("K1", "K2", "K3", "NK")

# Create an empty list to store the plots
plots_list <- list()

# Loop through each relatedness type
for (rel_type in relatedness_types) {

  # Create column names dynamically
  mean_col <- paste0("mean_", rel_type)
  psignif_col <- paste0("psignif_", rel_type, "_comb")
  prop_col <- paste0("mean_prop_", rel_type)

  # Create the plot title based on relatedness type
  if (rel_type == "K1") {
    x_label <- "Proportion of first-order kin mates"
  } else if (rel_type == "K2") {
    x_label <- "Proportion of second-order kin mates"
  } else if (rel_type == "K3") {
    x_label <- "Proportion of third-order kin mates"
  } else if (rel_type == "NK") {
    x_label <- "Proportion of non-kin mates"
  }

  # Generate the plot
  plot <- comb_BreederData_Sims_displim_nosexbias_PROPS %>%

    mutate(PairType = recode(PairType, First = "First pairings", Later = "Later pairings")) %>%
    mutate(datatype = factor(datatype, levels =c("Random mating with unconstrained dispersal", "Random mating with sex-biased limited dispersal", "Random mating with limited dispersal"))) %>%

    ### background distribution of relatedness of simulated pairs
    dplyr::filter(!datatype == "Observed") %>%
    dplyr::mutate(Sex = factor(Sex, levels = c("Female", "Male")),
                  PairType = factor(PairType, levels = c("First pairings", "Later pairings"))) %>%
    group_by(Sex, PairType) %>%
    ggplot(aes(x = !!sym(mean_col))) +
    geom_histogram(aes(alpha = !!sym(psignif_col), fill=datatype, color=datatype), bins = 30) +

    ## add vertical line of mean
    geom_vline(data = (dplyr::filter(comb_BreederData_Sims_displim_nosexbias_PROPS, datatype == "Observed") %>%
                         mutate(PairType = recode(PairType, First = "First pairings", Later = "Later pairings"))),
               color = "black", lwd=1.5, aes(xintercept = !!sym(mean_col))) +

    geom_segment(data = means_of_PROPS %>%
                   mutate(PairType = recode(PairType, First = "First pairings", Later = "Later pairings")),
                 aes(x = !!sym(prop_col), xend = !!sym(prop_col),
                     y = 400 * 1.5, yend = 400 * 1.2,
                     color = datatype),
                 arrow = arrow(length = unit(0.3, "cm")),
                 size = 1) +
    ## facet by sex and pairing type
    facet_wrap(vars(PairType, Sex), nrow=4, scales = "fixed") +

    ## labels and legends
    xlab(x_label) +
    ylab("Count") +

    # color, transparency and linetype aesthetics
    scale_color_viridis_d(option="inferno", begin = 0.2, end = 0.8) +
    scale_fill_viridis_d(option="inferno", begin = 0.2, end = 0.8, labels = label_wrap(19)) +
    scale_alpha_manual(values=c(0.5, 0)) +
    guides(alpha = "none", color = "none") +
    labs(fill = "Simulation type") +
    theme(legend.key.spacing.y = unit(0.1, "cm"))

  # Store the plot in the list
  plots_list[[rel_type]] <- plot
}

# Access individual plots (K2, K3, and NK commented out but can be run too)
K1 <- plots_list[["K1"]]
# K2 <- plots_list[["K2"]]
# K3 <- plots_list[["K3"]]
# NK <- plots_list[["NK"]]

# Display the plots (K2, K3, and NK commented out but can be run too)
K1
# K2
# K3
# NK

## Put two halves together for Figure 1 ===================================

library(patchwork)

p1 + K1 + plot_layout(guides = "collect", axes = "collect")



# 5. Fitness tradeoffs -----------------------------------------------------------------

FitnessData <- read_csv("./data/FitnessData.csv")

## Run models to asses the relationship between distance moved for first parings and weighted kinship coefficient on a breeders lifetime reproductive success. Breeding lifespan is included as a fixed effect to account for variation in lifespan. Natal Year is included as a random effect to account for birth cohort effects.

modelM <- glmer.nb(data=filter(FitnessData, Sex=="Male"), FitnessTotal ~ scale(DistFirstPairing) + scale(WtCoeff) + scale(BreedingLifespan) + (1|NatalYear))
summary(modelM) # Part of Table 1

modelM_sim_res <- simulateResiduals(modelM, n = 1000)
plotQQunif(modelM_sim_res)

modelF <- glmer.nb(data=filter(FitnessData, Sex=="Female"), FitnessTotal ~ scale(DistFirstPairing) + scale(WtCoeff) + scale(BreedingLifespan) + (1|NatalYear))
summary(modelF) # Part of Table 1

modelF_sim_res <- simulateResiduals(modelF, n = 1000)
plotQQunif(modelF_sim_res)

## For breeders that switched at least once in their lifetimes, run models to asses the relationship between total distance moved for later parings and weighted kinship coefficient on a breeders lifetime reproductive success. Breeding lifespan is included as a fixed effect to account for variation in lifespan. Natal Year is included as a random effect to account for birth cohort effects.

modelM_sub <- glmer.nb(data=filter(FitnessData, Sex=="Male" & MateSwitch==1), FitnessTotal ~ scale(TotalDistLaterPairing) + scale(WtCoeff) + scale(BreedingLifespan) + (1|NatalYear))
summary(modelM_sub) # Part of Table 2

modelM_sub_sim_res <- simulateResiduals(modelM_sub, n = 1000)
plotQQunif(modelM_sub_sim_res)

#modelF_sub <- glmer.nb(data=filter(FitnessData, Sex=="Female" & MateSwitch==1), FitnessTotal ~ scale(TotalDistLaterPairing) + scale(WtCoeff) + scale(BreedingLifespan) + (1|NatalYear))

## glmer.nb does not converge, use zero-inflated poisson instead

modelF_sub <- glmmTMB(
  FitnessTotal ~ scale(TotalDistLaterPairing) + scale(WtCoeff) + scale(BreedingLifespan) + (1|NatalYear),
  ziformula = ~ scale(BreedingLifespan),
  family = poisson(link = "log"),
  data = dplyr::filter(FitnessData, MateSwitch==1 & Sex=="Female")
)

summary(modelF_sub) # Part of Table 2

modelF_sub_sim_res <- simulateResiduals(modelF_sub, n = 1000)
plotQQunif(modelF_sub_sim_res)

## Make Fig. 2 =====================================

LRS_First_Dist <- ggplot(FitnessData, aes(x = DistFirstPairing, y = FitnessTotal, group = Sex)) +
                      geom_jitter(aes(shape=Sex, color=Sex), height=0.5, width=0.20, alpha=0.3) +
                      geom_smooth(method = "glm.nb", aes(color=Sex, fill=Sex, linetype=Sex), alpha=0.2,
                                  show.legend = c(color=TRUE, fill=TRUE, linetype=FALSE))  +
  					scale_color_manual(values = c("orange", "blue3")) +
                      scale_fill_manual(values = c("orange", "blue3")) +
                      scale_linetype_manual(values = c("dashed", "dashed"),drop = FALSE) +
                      xlab("Distance moved \n for first pairing") +
                      xlim(-0.20, 15) +
                      ylab("") +
  					guides(color = guide_legend(override.aes = list(shape = c(17, 16), size=2, linetype="solid"))) +
                   scale_shape(guide = 'none')

LRS_First_Coeff <- ggplot(FitnessData, aes(x = WtCoeff, y = FitnessTotal, group = Sex)) +
  				  	geom_jitter(aes(shape=Sex, color=Sex), height=0.5, width=0.005, alpha=0.3) +
   					geom_smooth(method = "glm.nb", aes(color=Sex, fill=Sex, linetype=Sex), alpha=0.2,
               				 show.legend = c(color=TRUE, fill=TRUE, linetype=FALSE)) +
  					scale_color_manual(values = c("orange", "blue3")) +
  					scale_fill_manual(values = c("orange", "blue3")) +
  					scale_linetype_manual(values = c("dashed", "solid"),drop = FALSE) +
  					xlab("Weighted kinship \n coefficient to mates") +
  					xlim(-0.005, (0.25+0.005)) +
  					ylab("") +
  					guides(color = guide_legend(override.aes = list(shape = c(17, 16), size=2, linetype="solid"))) +
  					scale_shape(guide = 'none')

## Make model predictions for modelM_sub and modelF_sub since I'm using negative binomial for the former and zero-inflated poisson for the latter

pred_WtCoeff_F <- ggpredict(modelF_sub, terms = "WtCoeff", bias_correction = TRUE)
pred_WtCoeff_F$Sex <- "Female"
pred_WtCoeff_M <- ggpredict(modelM_sub, terms = "WtCoeff", bias_correction = TRUE)
pred_WtCoeff_M$Sex <- "Male"

pred_WtCoeff <- bind_rows(pred_WtCoeff_F, pred_WtCoeff_M)

pred_WtCoeff <- as.data.frame(pred_WtCoeff) %>%
				 dplyr::rename(WtCoeff = x,
    						   FitnessTotal = predicted)

pred_TotalDistLaterPairing_F <- ggpredict(modelF_sub, terms = "TotalDistLaterPairing")
pred_TotalDistLaterPairing_F$Sex <- "Female"
pred_TotalDistLaterPairing_M <- ggpredict(modelM_sub, terms = "TotalDistLaterPairing")
pred_TotalDistLaterPairing_M$Sex <- "Male"

pred_TotalDistLaterPairing <- bind_rows(pred_TotalDistLaterPairing_F, pred_TotalDistLaterPairing_M)

pred_TotalDistLaterPairing <- as.data.frame(pred_TotalDistLaterPairing) %>%
  						   dplyr::rename(TotalDistLaterPairing = x,
  						   				  FitnessTotal = predicted)

LRS_Later_TotalDistLaterPairing <- ggplot(filter(FitnessData, MateSwitch==1),
											 aes(x = TotalDistLaterPairing, y = FitnessTotal, group = Sex)) +
  									  geom_jitter(aes(shape=Sex, color=Sex), height=0.5, width=0.20, alpha=0.3) +
  									  geom_ribbon(data= pred_TotalDistLaterPairing,
  									  			   aes(ymin = conf.low, ymax = conf.high, fill = Sex),
  									  			   alpha = 0.2) +
  									  geom_line(data= pred_TotalDistLaterPairing,
  									  			 aes(color = Sex, linetype=Sex),
  									  			 size = 1) +
  									  scale_color_manual(values = c("orange", "blue3")) +
  									  scale_fill_manual(values = c("orange", "blue3")) +
  scale_linetype_manual(values = c("dashed", "dashed"),drop = FALSE) +
  xlab("Total distance moved \n for later pairings") +
  xlim(-0.20, 15) +
  ylab("") +
  guides(color = guide_legend(override.aes = list(shape = c(17, 16), size=2, linetype="solid"))) +
  scale_shape(guide = 'none')

LRS_Later_Coeff <- ggplot(filter(FitnessData, MateSwitch==1),
							  aes(x = WtCoeff, y = FitnessTotal, group = Sex)) +
  geom_jitter(aes(shape=Sex, color=Sex), height=0.5, width=0.005, alpha=0.3) +
  geom_ribbon(data= pred_WtCoeff, aes(ymin = conf.low, ymax = conf.high, fill = Sex), alpha = 0.2) +
  geom_line(data= pred_WtCoeff, aes(color = Sex, linetype=Sex), size = 1) +
  scale_color_manual(values = c("orange", "blue3")) +
  scale_fill_manual(values = c("orange", "blue3")) +
  scale_linetype_manual(values = c("dashed", "solid"), drop = FALSE) +
  xlab("Weighted kinship \n coefficient to mates") +
  xlim(-0.005, (0.25+0.005)) +
  ylab("") +
  guides(color = guide_legend(override.aes = list(shape = c(17, 16), size=2, linetype="solid"))) +
  scale_shape(guide = 'none')


## Put plots together for Fig. 2

((LRS_First_Dist + LRS_First_Coeff) / (LRS_Later_TotalDistLaterPairing + LRS_Later_Coeff)) + plot_layout(axis_titles = "collect_y", guides = "collect") + plot_annotation(tag_levels = 'A')

# 6. Survival --------------------------

SurvivalData <- read_csv("./data/SurvivalData.csv")

## Fit Cox Proportional Hazards models to assess the relationship between distance moved for first pairings (natal dispersal distance) and loss of breeding territory during later pairings (site fidelity) on breeder survival

## Natal dispersal distance.

SurvResults_split <- list()
exit_threshold <- 4 # fit stratified Cox proportional hazards models with two time intervals: ≤ 4 years and > 4 years from the start of breeding to account for a time-varying effect of distance moved such that the proportional hazard assumption was met (see SI for details).

  for (sex in c("Male", "Female")) {
    df_sex <- subset(SurvivalData, Sex == sex)
    df_sex <- df_sex %>% mutate(tgroup = ifelse(ExitAge <= exit_threshold, 1, 2))
    model <- coxph(Surv(EntryAge, ExitAge, Dead) ~ Dist:strata(tgroup) + cluster(NatalYear), data = df_sex)
    test.hazards <- cox.zph(model)
    print(test.hazards)  # CHECK PROPORTIONAL HAZARD ASSUMPTION
    name <- paste("mSurvival", sex, "ExitThreshold", exit_threshold, sep = "_")
    SurvResults_split[[name]] <- model
  }

SurvResults_split

## Site fidelity

SurvResults_sitefidelity <- list()

for (sex in c("Male", "Female")) {
  df_sex <- subset(SurvivalData, Sex==sex & EntryAge>0 & numMateSwitch>0)
  model <- survival::coxph(survival::Surv(EntryAge, ExitAge, Dead) ~ as.factor(TerrLost) + cluster(NatalYear), data = df_sex)
  test.hazards <- (survival::cox.zph(model))
  print(test.hazards) #CHECK PROPORTIONAL HAZARD ASSUMPTION
  print(paste("Sample size:", length(unique(df_sex$FocalInd))))
  SurvResults_sitefidelity[[paste("mSurvival", sex, sep = "_")]]  <- model
}

SurvResults_sitefidelity

# 7. Divorce --------------------------

## Fit GLMMs to assess the relationship between kinship coefficient to mate and likelihood of divorce. Separate models by sex and pairing type.

summary(glmer(Divorce ~ scale(Coeff) + (1|Year), data=filter(BreederData, Sex=="Male" & PairType=="First"), family="binomial"))

summary(glmer(Divorce ~ scale(Coeff) + (1|Year) + (1|FocalInd), data=filter(BreederData, Sex=="Male" & PairType=="Later"), family="binomial"))

summary(glmer(Divorce ~ scale(Coeff) + (1|Year), data=filter(BreederData, Sex=="Female" & PairType=="First"), family="binomial"))

summary(glmer(Divorce ~ scale(Coeff) + (1|Year) + (1|FocalInd), data=filter(BreederData, Sex=="Female" & PairType=="Later"), family="binomial"))

## Fit GLMMs to assess the relationship between the mate being a first order relative and likelihood of divorce. Separate models by sex and pairing type.

summary(glmer(Divorce ~ as.factor(FirstOrderRel) + (1|Year), data=filter(BreederData, Sex=="Male" & PairType=="First"), family="binomial"))

summary(glmer(Divorce ~ as.factor(FirstOrderRel) + (1|Year) + (1|FocalInd), data=filter(BreederData, Sex=="Male" & PairType=="Later"), family="binomial"))

summary(glmer(Divorce ~ as.factor(FirstOrderRel) + (1|Year), data=filter(BreederData, Sex=="Female" & PairType=="First"), family="binomial"))

summary(glmer(Divorce ~ as.factor(FirstOrderRel) + (1|Year) + (1|FocalInd), data=filter(BreederData, Sex=="Female" & PairType=="Later"), family="binomial"))
