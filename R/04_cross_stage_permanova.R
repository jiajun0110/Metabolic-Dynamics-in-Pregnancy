## Cross-stage PERMANOVA and Spearman correlations
## Input: early and late metabolite-change matrices, with matched rows.

library(vegan)

cross_stage_permanova <- function(early, late, seed) {
  early <- as.data.frame(early)
  dist_late <- dist(as.matrix(late), method = "euclidean")
  set.seed(seed)
  result <- lapply(names(early), function(met) {
    predictor <- data.frame(change = early[[met]])
    fit <- adonis2(dist_late ~ change, data = predictor, permutations = 999)
    data.frame(metabolite = met, R2 = fit$R2[1], P = fit$`Pr(>F)`[1])
  })
  result <- do.call(rbind, result)
  result$FDR <- p.adjust(result$P, "BH")
  result
}

cross_stage_correlations <- function(early, late) {
  early <- as.data.frame(early)
  late <- as.data.frame(late)
  pairs <- expand.grid(early = names(early), late = names(late), stringsAsFactors = FALSE)
  result <- lapply(seq_len(nrow(pairs)), function(i) {
    test <- cor.test(early[[pairs$early[i]]], late[[pairs$late[i]]],
                     method = "spearman", exact = FALSE)
    data.frame(pairs[i, ], rho = unname(test$estimate), P = test$p.value)
  })
  result <- do.call(rbind, result)
  result$FDR <- p.adjust(result$P, "BH")
  result
}
