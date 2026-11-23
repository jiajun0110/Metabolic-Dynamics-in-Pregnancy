## Fetal growth rates, metabolite associations and tertile comparisons.
## Inputs: analysis-ready ultrasound (id, gestational_week, efw) and maternal data.
## Supply metabolite/covariate columns and group as a T1/T2/T3 factor, with T1 first.

fit_fetal_growth <- function(ultrasound) {
  model <- lmerTest::lmer(efw ~ gestational_week + (gestational_week | id),
                          data = ultrasound)
  slopes <- lme4::ranef(model)$id
  rates <- data.frame(id = rownames(slopes),
    weekly_growth = unname(lme4::fixef(model)["gestational_week"]) +
      slopes$gestational_week)
  list(model = model, rates = rates)
}

fetal_growth_correlations <- function(data, metabolites,
    covariates = c("age", "parity", "gravidity", "prebmi"), discovery = TRUE) {
  results <- lapply(metabolites, function(metabolite) {
    test <- ppcor::pcor.test(data$weekly_growth, data[[metabolite]],
                            data[, covariates, drop = FALSE], method = "spearman")
    rho <- unname(test$estimate)
    ci <- tanh(atanh(rho) + c(-1, 1) * 1.96 /
                 sqrt(nrow(data) - 3 - length(covariates)))
    data.frame(metabolite = metabolite, n = nrow(data), rho = rho,
               ci_lower = ci[1], ci_upper = ci[2], p_value = test$p.value)
  })
  results <- do.call(rbind, results)
  if (discovery) results$fdr <- stats::p.adjust(results$p_value, method = "BH")
  results
}

fit_fetal_growth_tertiles <- function(ultrasound, growth_rates) {
  visits <- droplevels(ultrasound[ultrasound$group %in% c("T1", "T3"), ])
  rates <- droplevels(growth_rates[growth_rates$group %in% c("T1", "T3"), ])
  model <- lmerTest::lmer(
    efw ~ gestational_week * group + (gestational_week | id), data = visits)
  test <- stats::wilcox.test(weekly_growth ~ group, data = rates,
                             alternative = "two.sided", exact = FALSE)
  interaction <- coef(summary(model))["gestational_week:groupT3", , drop = FALSE]
  list(model = model, interaction = interaction, wilcoxon = test)
}

fit_fetal_fgr <- function(data) {
  model <- stats::glm(fgr ~ group + age + parity + gravidity + prebmi,
                      data = data, family = stats::binomial())
  coefficients <- coef(summary(model))[c("groupT2", "groupT3"), , drop = FALSE]
  ci <- suppressMessages(stats::confint(model, parm = c("groupT2", "groupT3")))
  results <- data.frame(group = c("T2", "T3"),
    odds_ratio = exp(coefficients[, "Estimate"]),
    ci_lower = exp(ci[, 1]), ci_upper = exp(ci[, 2]),
    p_value = coefficients[, "Pr(>|z|)"])
  list(model = model, results = results)
}
