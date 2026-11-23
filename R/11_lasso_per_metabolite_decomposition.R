## Per-metabolite adjusted-R2 decomposition.
## Inputs: row-aligned y and factor_list named genetics, microbiome, clinical, dietary.
## Genetics contains GWAS-selected variants; dietary variables are prepared factors.
## Predictor names are unique and syntactically valid. R2 values are fractions.

lasso_select <- function(y, X, seed) {
  data <- stats::na.omit(data.frame(y = y, X))
  X <- droplevels(data[, -1, drop = FALSE])
  X <- X[, vapply(X, function(x) length(unique(x)) > 1, logical(1)), drop = FALSE]
  if (!ncol(X)) return(character())
  matrix <- stats::model.matrix(~ ., data = X)
  variable_index <- attr(matrix, "assign")[-1]
  matrix <- matrix[, -1, drop = FALSE]
  if (ncol(matrix) < 2) return(character())
  set.seed(seed)
  fit <- glmnet::cv.glmnet(matrix, data$y, alpha = 1, nfolds = 10)
  beta <- as.matrix(stats::coef(fit, s = "lambda.min"))
  columns <- setdiff(rownames(beta)[beta[, 1] != 0], "(Intercept)")
  unique(names(X)[variable_index[match(columns, colnames(matrix))]])
}

decompose_one <- function(y, factor_list, seed) {
  # Select each non-genetic domain independently.
  selected <- lapply(names(factor_list), function(group) {
    if (group == "genetics") names(factor_list[[group]]) else
      lasso_select(y, factor_list[[group]], seed)
  })
  names(selected) <- names(factor_list)
  X <- do.call(cbind, unname(Map(function(data, variables) data[, variables, drop = FALSE],
                                 factor_list, selected)))
  data <- stats::na.omit(data.frame(y = y, X))
  full <- stats::lm(y ~ ., data = data)
  full_r2 <- summary(full)$adj.r.squared

  # Full and reduced models use the same complete observations.
  category_r2 <- vapply(selected, function(variables) {
    reduced <- stats::lm(y ~ ., data = data[, setdiff(names(data), variables), drop = FALSE])
    full_r2 - summary(reduced)$adj.r.squared
  }, numeric(1))
  list(model = full, selected = selected, adjusted_R2 = full_r2,
       category_R2 = category_r2)
}
