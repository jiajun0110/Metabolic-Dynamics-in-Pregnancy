## Backward-selection PERMANOVA across factor categories.
## Inputs: complete, row-aligned response matrix y and named factor_list.
## Factor columns have unique, syntactically valid names. R2 values are fractions.

backward_permanova <- function(y, factor_list, seed, nperm = 999) {
  X <- do.call(cbind, unname(factor_list))
  category <- setNames(rep(names(factor_list), vapply(factor_list, ncol, integer(1))), names(X))
  distance <- stats::dist(y, method = "euclidean")

  # Screen individual factors.
  p <- vapply(names(X), function(variable) {
    set.seed(seed)
    fit <- vegan::adonis2(stats::reformulate(variable, response = "distance"),
                          data = X, permutations = nperm)
    fit$`Pr(>F)`[1]
  }, numeric(1))
  selected <- names(p)[!is.na(p) & p < 0.05]
  category_r2 <- setNames(rep(0, length(factor_list)), names(factor_list))
  model <- NULL

  # Remove the least significant marginal term until all retained terms pass.
  while (length(selected)) {
    set.seed(seed)
    model <- vegan::adonis2(stats::reformulate(selected, response = "distance"),
                            data = X, permutations = nperm, by = "margin")
    marginal_p <- model[selected, "Pr(>F)"]
    if (all(marginal_p < 0.05)) break
    selected <- selected[-which.max(marginal_p)]
  }
  if (!length(selected)) return(list(model = NULL, selected = selected,
                                     category_R2 = category_r2, total_R2 = 0))

  for (group in names(factor_list)) {
    terms <- selected[category[selected] == group]
    category_r2[group] <- sum(model[terms, "R2"])
  }
  list(model = model, selected = selected, category_R2 = category_r2,
       total_R2 = 1 - model["Residual", "R2"])
}
