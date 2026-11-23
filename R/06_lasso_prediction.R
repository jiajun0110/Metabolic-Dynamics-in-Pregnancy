## LASSO and XGBoost prediction of continuous outcomes
## Inputs: prepared numeric training/validation tables and column names.
## Supply XGBoost params, nrounds and random seeds.

run_lasso <- function(X, y, seed) {
  set.seed(seed)
  foldid <- sample(rep(seq_len(10), length.out = nrow(X)))
  fit <- glmnet::cv.glmnet(as.matrix(X), y, alpha = 1, family = "gaussian",
                          nfolds = 10, foldid = foldid)
  coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
  setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
}

bootstrap_spearman <- function(y, pred_crf, pred_full, seed, R = 2000) {
  set.seed(seed)
  rho <- c(stats::cor(y, pred_crf, method = "spearman"),
           stats::cor(y, pred_full, method = "spearman"))
  indices <- matrix(sample.int(length(y), length(y) * R, replace = TRUE), nrow = length(y))
  draws <- t(vapply(seq_len(R), function(b) {
    i <- indices[, b]
    c(stats::cor(y[i], pred_crf[i], method = "spearman"),
      stats::cor(y[i], pred_full[i], method = "spearman"))
  }, numeric(2)))
  colnames(draws) <- c("CRF", "CRF_Metabolites")
  diffs <- draws[, 2] - draws[, 1]
  delta <- rho[2] - rho[1]
  p <- (1 + sum(abs(diffs - delta) >= abs(delta))) / (length(diffs) + 1)
  list(rho_crf = rho[1], rho_full = rho[2], delta_rho = delta,
       rho_ci = apply(draws, 2, stats::quantile, probs = c(0.025, 0.975)),
       delta_ci = stats::quantile(diffs, c(0.025, 0.975)), p_value = p, bootstrap = draws)
}

# For sampling-time sensitivity, supply restricted data and the main selection.
compare_models <- function(train, valid, outcome, crf_vars, metabolite_vars,
                           params, nrounds, seed, selected_metabolites = NULL,
                           bootstrap_seed = seed, R = 2000) {
  set.seed(seed)
  predictors <- c(crf_vars, metabolite_vars)
  columns <- c(outcome, predictors)
  train <- train[stats::complete.cases(train[, columns]), , drop = FALSE]
  valid <- valid[stats::complete.cases(valid[, columns]), , drop = FALSE]
  center <- colMeans(train[, predictors, drop = FALSE])
  spread <- vapply(train[, predictors, drop = FALSE], stats::sd, numeric(1))
  spread[spread == 0] <- 1
  X_train <- scale(as.matrix(train[, predictors, drop = FALSE]), center, spread)
  X_valid <- scale(as.matrix(valid[, predictors, drop = FALSE]), center, spread)
  if (is.null(selected_metabolites)) {
    selected_metabolites <- intersect(run_lasso(X_train, train[[outcome]], seed), metabolite_vars)
  }
  full_vars <- c(crf_vars, selected_metabolites)
  set.seed(seed)
  m1 <- xgboost::xgb.train(params = params, nrounds = nrounds, verbose = 0,
    data = xgboost::xgb.DMatrix(X_train[, crf_vars, drop = FALSE], label = train[[outcome]]))
  set.seed(seed)
  m2 <- xgboost::xgb.train(params = params, nrounds = nrounds, verbose = 0,
    data = xgboost::xgb.DMatrix(X_train[, full_vars, drop = FALSE], label = train[[outcome]]))
  pred1 <- as.numeric(stats::predict(m1, X_valid[, crf_vars, drop = FALSE]))
  pred2 <- as.numeric(stats::predict(m2, X_valid[, full_vars, drop = FALSE]))
  list(model_crf = m1, model_full = m2, selected_metabolites = selected_metabolites,
       center = center, scale = spread, train_X = X_train[, full_vars, drop = FALSE],
       predictions = data.frame(observed = valid[[outcome]], pred_crf = pred1, pred_full = pred2),
       performance = bootstrap_spearman(valid[[outcome]], pred1, pred2, bootstrap_seed, R))
}

signed_shap <- function(model, X) {
  X <- as.matrix(X)
  values <- stats::predict(model, X, predcontrib = TRUE, approxcontrib = FALSE)
  values <- values[, seq_len(ncol(X)), drop = FALSE]
  direction <- vapply(seq_len(ncol(X)), function(j)
    sign(suppressWarnings(stats::cor(X[, j], values[, j], method = "spearman"))), numeric(1))
  direction[is.na(direction)] <- 0
  importance <- colMeans(abs(values))
  data.frame(feature = colnames(X), mean_abs_shap = importance,
             signed_mean_abs_shap = importance * direction, row.names = NULL)
}
