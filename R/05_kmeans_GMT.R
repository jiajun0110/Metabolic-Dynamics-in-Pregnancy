## K-means clustering for gestational metabolic trajectory (GMT) types
## Input: numeric matrix of T1-to-T2 metabolite changes.

library(NbClust)
library(fpc)

fit_gmt <- function(dynamics_matrix, seed) {
  X <- as.matrix(dynamics_matrix)
  bounds <- apply(X, 2, quantile, probs = c(0.01, 0.99))
  for (j in seq_len(ncol(X))) X[, j] <- pmin(pmax(X[, j], bounds[1, j]), bounds[2, j])
  set.seed(seed)
  nb <- NbClust(X, min.nc = 2, max.nc = 6, method = "kmeans", index = "all")
  set.seed(seed)
  km <- kmeans(X, centers = 2, nstart = 25, iter.max = 100)
  stability <- clusterboot(X, B = 1000, bootmethod = "boot",
                            clustermethod = kmeansCBI, krange = 2, seed = seed)
  list(model = km, centers = km$centers, bounds = bounds, selection = nb, stability = stability)
}

assign_gmt <- function(dynamics_matrix, fit) {
  X <- as.matrix(dynamics_matrix)[, colnames(fit$centers), drop = FALSE]
  for (j in seq_len(ncol(X))) X[, j] <- pmin(pmax(X[, j], fit$bounds[1, j]), fit$bounds[2, j])
  distance <- do.call(cbind, lapply(seq_len(nrow(fit$centers)), function(k)
    rowSums(sweep(X, 2, fit$centers[k, ], "-")^2)))
  max.col(-distance, ties.method = "first")
}

compare_gmt_changes <- function(changes, group, log2fc) {
  score <- rowMeans(scale(as.matrix(changes)))
  p <- vapply(as.data.frame(log2fc), function(x) wilcox.test(x ~ group)$p.value, numeric(1))
  list(score = score, overall_test = wilcox.test(score ~ group),
       metabolites = data.frame(metabolite = names(p), P = p, FDR = p.adjust(p, "BH")))
}
