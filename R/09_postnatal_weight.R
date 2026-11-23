## Postnatal weight-for-age z-score and weight growth models.
## Inputs: analysis-ready follow-up with month, waz/weight, metabolite and covariates.
## WAZ uses child_id; weight uses id and a prepared T1/T2/T3 group factor, T1 first.

fit_postnatal_waz <- function(data, metabolite_sd) {
  model <- lmerTest::lmer(
    waz ~ month + metabolite + age + prebmi + parity + gravidity + sex + bw_z +
      (month | child_id), data = data,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
  coefficient <- coef(summary(model))["metabolite", ]
  beta <- unname(coefficient["Estimate"])
  ci <- (beta + c(-1, 1) * stats::qt(0.975, coefficient["df"]) *
           coefficient["Std. Error"]) * metabolite_sd
  results <- data.frame(beta_raw = beta, beta_per_sd = beta * metabolite_sd,
                        ci_lower = ci[1], ci_upper = ci[2],
                        p_value = unname(coefficient["Pr(>|t|)"]))
  list(model = model, results = results)
}

fit_postnatal_weight <- function(data) {
  data <- data[data$month >= 1 & data$month <= 36, ]
  growth_model <- lmerTest::lmer(weight ~ month + (month | id), data = data)
  slopes <- lme4::ranef(growth_model)$id
  rates <- data.frame(id = rownames(slopes),
    monthly_gain = unname(lme4::fixef(growth_model)["month"]) + slopes$month,
    group = data$group[match(rownames(slopes), data$id)])
  extreme_rates <- droplevels(rates[rates$group %in% c("T1", "T3"), ])
  test <- stats::wilcox.test(monthly_gain ~ group, data = extreme_rates,
                             alternative = "two.sided", exact = FALSE)
  visits <- droplevels(data[data$group %in% c("T1", "T3"), ])
  trajectory_model <- lmerTest::lmer(weight ~ month * group + (month | id),
                                     data = visits)
  interaction <- coef(summary(trajectory_model))["month:groupT3", , drop = FALSE]
  list(growth_model = growth_model, rates = rates, trajectory_model = trajectory_model,
       interaction = interaction, wilcoxon = test)
}
