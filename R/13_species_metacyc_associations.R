## Species-metabolite and species-stratified MetaCyc associations.
## Input data contains metabolite changes, prepared species/pathway abundances,
## and age, parity, gravidity, prebmi. pairs specifies the selected column pairs.

species_associations <- function(data, pairs) {
  # pairs has metabolite and species columns; BH spans all supplied pairs.
  result <- lapply(seq_len(nrow(pairs)), function(i) {
    d <- data[, c(pairs$metabolite[i], pairs$species[i], "age", "parity", "gravidity", "prebmi")]
    names(d)[1:2] <- c("metabolite", "species")
    fit <- stats::lm(metabolite ~ species + age + parity + gravidity + prebmi, data = d)
    coefficient <- summary(fit)$coefficients["species", ]
    data.frame(pairs[i, ], beta = coefficient["Estimate"],
               SE = coefficient["Std. Error"], P = coefficient["Pr(>|t|)"], row.names = NULL)
  })
  result <- do.call(rbind, result)
  result$FDR <- stats::p.adjust(result$P, method = "BH")
  result
}

metacyc_associations <- function(data, pairs) {
  # pairs has metabolite and pathway columns for validated species-metabolite pairs.
  result <- lapply(seq_len(nrow(pairs)), function(i) {
    test <- stats::cor.test(data[[pairs$metabolite[i]]], data[[pairs$pathway[i]]],
                            method = "spearman", exact = FALSE)
    data.frame(pairs[i, ], rho = unname(test$estimate), P = test$p.value)
  })
  result <- do.call(rbind, result)
  result$FDR <- stats::p.adjust(result$P, method = "BH")
  result
}
