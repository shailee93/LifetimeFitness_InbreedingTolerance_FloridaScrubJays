# Potential mate kinship coeffiicent with distance analysis

PotentialMateCoeffFn <- function(df) {

  require(lme4)
  require(jtools)
  require(ggpattern)
  #require(sjPlot)

  # Define the vectors
  PairType <- c("First", "Later")
  Sex <- c("Male", "Female")

  # Create combinations of phases and sexes
  combinations <- expand.grid(PairType = PairType, Sex = Sex)

  # Run models for each combination using lapply
  models <- lapply(1:nrow(combinations), function(i) {
    glmer(
      formula = CloseRel ~ scale(Dist) + (1|Year) + (1|FocalInd),
      data = dplyr::filter(df, PairType == combinations$PairType[i] & Sex == combinations$Sex[i]),
      family = "binomial"
    )
  })

  # Summarize the results
  summaries <- lapply(models, jtools::summ, confint=TRUE, digits=3, exp=TRUE)

  return(list(summaries))
}
