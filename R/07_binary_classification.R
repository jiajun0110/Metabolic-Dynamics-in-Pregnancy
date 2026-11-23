## Logistic classification and GMT clinical associations
## Inputs: prepared tables, binary outcomes coded 0/1, and predictor column names.
## Include disease-history variables in crf_vars for disease-specific comparisons.

bootstrap_auc <- function(y, pred_crf, pred_full, seed, R = 2000) {
  set.seed(seed)
  r1 <- pROC::roc(y, pred_crf, levels = c(0, 1), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(y, pred_full, levels = c(0, 1), direction = "<", quiet = TRUE)
  test <- pROC::roc.test(r1, r2, paired = TRUE, method = "bootstrap",
                        boot.n = R, boot.stratified = TRUE, alternative = "two.sided")
  control <- which(y == 0); cases <- which(y == 1)
  set.seed(seed)
  draws <- t(vapply(seq_len(R), function(b) {
    i <- c(control[sample.int(length(control), length(control), replace = TRUE)],
           cases[sample.int(length(cases), length(cases), replace = TRUE)])
    c(as.numeric(pROC::auc(y[i], pred_crf[i], levels = c(0, 1), direction = "<", quiet = TRUE)),
      as.numeric(pROC::auc(y[i], pred_full[i], levels = c(0, 1), direction = "<", quiet = TRUE)))
  }, numeric(2)))
  colnames(draws) <- c("CRF", "CRF_Metabolites")
  aucs <- c(as.numeric(pROC::auc(r1)), as.numeric(pROC::auc(r2)))
  delta <- aucs[2] - aucs[1]
  D <- as.numeric(test$statistic)
  se <- if (is.finite(D) && D != 0) abs(delta / D) else stats::sd(draws[, 2] - draws[, 1])
  p <- if (delta == 0 && se == 0) 1 else as.numeric(test$p.value)
  list(auc_crf = aucs[1], auc_full = aucs[2], delta_auc = delta,
       auc_ci = apply(draws, 2, stats::quantile, probs = c(0.025, 0.975)),
       delta_ci = delta + c(-1, 1) * stats::qnorm(0.975) * se, p_value = p, bootstrap = draws)
}

compare_binary <- function(data, outcome, crf_vars, metabolite_vars, seed, R = 2000) {
  set.seed(seed)
  columns <- c(outcome, crf_vars, metabolite_vars)
  d <- data[stats::complete.cases(data[, columns]), , drop = FALSE]
  m1 <- stats::glm(stats::reformulate(crf_vars, outcome), data = d, family = stats::binomial())
  m2 <- stats::glm(stats::reformulate(c(crf_vars, metabolite_vars), outcome),
                   data = d, family = stats::binomial())
  p1 <- as.numeric(stats::fitted(m1)); p2 <- as.numeric(stats::fitted(m2))
  list(model_crf = m1, model_full = m2,
       predictions = data.frame(observed = d[[outcome]], pred_crf = p1, pred_full = p2),
       performance = bootstrap_auc(d[[outcome]], p1, p2, seed, R))
}

gmt_clinical_associations <- function(data, group_var, continuous, binary, crf_vars) {
  continuous_results <- lapply(continuous, function(y) {
    test <- stats::wilcox.test(stats::reformulate(group_var, y), data = data)
    data.frame(outcome = y, p_value = test$p.value)
  })
  binary_models <- lapply(binary, function(y)
    stats::glm(stats::reformulate(c(group_var, crf_vars), y), data = data, family = stats::binomial()))
  list(continuous = do.call(rbind, continuous_results),
       binary = stats::setNames(binary_models, binary))
}
