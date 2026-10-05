# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

ObsSimsFn <- function(ObsDF, SimDF, cVars, cSexes, cPairType, StartYear, EndYear, EmpPValueType) {

  # Prepare the observed data
  ObsDF <- ObsDF %>%
    dplyr::select(c(Year, Sex, PairType, cVars)) %>%
    dplyr::mutate(DataType = "Obs") %>%
    tidyr::pivot_longer(cols = cVars, names_to = "Var", values_to = "value")

  # Prepare the simulated data
  SimDF <- SimDF %>%
    dplyr::select(c(Year, Sex, PairType, Rep, cVars)) %>%
    dplyr::mutate(DataType = "Sim") %>%
    tidyr::pivot_longer(cols = cVars, names_to = "Var", values_to = "value")

  # Group and summarize
    ObsSimMeansDF <- bind_rows(ObsDF, SimDF) %>%
      group_by(DataType, PairType, Sex, Rep, Var) %>%
      summarize(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
      ungroup()

l <- list()
values <- list()
Sex <- list()
PairType <- list()
VarL <- list()
counter <- 0

for (a in cVars) {
  print(a)
  ObsSimMean_a <- filter(ObsSimMeansDF, Var == a)

  for (d in cSexes) {
    print(d)
    ObsSimMean_d <- dplyr::filter(ObsSimMean_a, Sex == d)
    ObsSimMean_g <- ObsSimMean_d

      for (f in cPairType) {
        print(f)
        ObsSimMean_f <- dplyr::filter(ObsSimMean_g, PairType == f)

        if (!nrow(dplyr::filter(ObsSimMean_f, DataType == "Obs")) == 0 &
            !nrow(dplyr::filter(ObsSimMean_f, DataType == "Sim")) == 0) {
          ObsMeanVar <- dplyr::filter(ObsSimMean_f, DataType == "Obs")$mean
          SimMeanVar <- dplyr::filter(ObsSimMean_f, DataType == "Sim")$mean

          if (EmpPValueType == "more") {
            S <- sum(SimMeanVar >= ObsMeanVar)
          }

          if (EmpPValueType == "less") {
            S <- sum(SimMeanVar <= ObsMeanVar)
          }

          counter <- counter + 1
          N <- nrow(dplyr::filter(ObsSimMean_f, DataType == "Sim"))
          p <- (S+1) / (N+1)
        } else {
          p <- NA
          ObsMeanVar <- NA
        }

        l[[counter]] <- p
        values[[counter]] <- ObsMeanVar
        Sex[[counter]] <- d
        VarL[[counter]] <- a
        PairType[[counter]] <- f
      }
    }
  }



outputDF <- data.frame(Var = unlist(VarL),
                       ObsMean = unlist(values),
                       p = unlist(l),
                       Sex = unlist(Sex),
                       PairType = unlist(PairType))

finalDF <- merge(ObsSimMeansDF, outputDF, by = c("Sex", "Var", "PairType"), all.x = TRUE, all.y = TRUE)

return(finalDF)
}
