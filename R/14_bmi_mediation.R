## BMI-metabolite associations and mediation of the BMI-OGTT2h association.
## Input data contains prebmi, age, parity, gravidity and metabolite changes.
## Mediation uses generic columns metabolite (the selected change) and ogtt2h.

bmi_associations <- function(data, metabolites) {
  result <- lapply(metabolites, function(variable) {
    d <- data[, c(variable, "prebmi", "age", "parity", "gravidity")]
    names(d)[1] <- "metabolite"
    fit <- stats::lm(metabolite ~ prebmi + age + parity + gravidity, data = d)
    coefficient <- summary(fit)$coefficients["prebmi", ]
    data.frame(metabolite = variable, beta = coefficient["Estimate"],
               SE = coefficient["Std. Error"], P = coefficient["Pr(>|t|)"], row.names = NULL)
  })
  result <- do.call(rbind, result)
  result$FDR <- stats::p.adjust(result$P, method = "BH")
  result
}

run_mediation <- function(data, seed, B = 2000) {
  d <- stats::na.omit(data[, c("metabolite", "prebmi", "ogtt2h", "age", "parity", "gravidity")])
  mediator_model <- stats::lm(metabolite ~ prebmi + age + parity + gravidity, data = d)
  outcome_model <- stats::lm(ogtt2h ~ prebmi + metabolite + age + parity + gravidity, data = d)
  set.seed(seed)
  mediation::mediate(mediator_model, outcome_model, treat = "prebmi",
                      mediator = "metabolite", boot = TRUE, sims = B)
}
