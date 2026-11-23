## I² heterogeneity for stage-exclusive metabolites
##
## Input: paired_fc — data.frame [id, metabolite, fc_T1_to_T2, fc_T2_to_T3]
##        (log2 fold changes for each individual across both transitions)

library(dplyr)

calc_I2 <- function(paired_fc) {
  paired_fc %>%
    group_by(metabolite) %>%
    summarise(
      mean_early = mean(fc_T1_to_T2, na.rm = TRUE),
      se_early   = sd(fc_T1_to_T2, na.rm = TRUE) / sqrt(sum(!is.na(fc_T1_to_T2))),
      mean_late  = mean(fc_T2_to_T3, na.rm = TRUE),
      se_late    = sd(fc_T2_to_T3, na.rm = TRUE) / sqrt(sum(!is.na(fc_T2_to_T3))),
      .groups = "drop"
    ) %>%
    mutate(
      w1 = 1 / se_early^2, w2 = 1 / se_late^2,
      summary_effect = (w1 * mean_early + w2 * mean_late) / (w1 + w2),
      Q  = w1 * (mean_early - summary_effect)^2 + w2 * (mean_late - summary_effect)^2,
      I2 = pmax(0, (Q - 1) / Q) * 100,
      p_heterogeneity = pchisq(Q, df = 1, lower.tail = FALSE)
    ) %>%
    mutate(FDR = p.adjust(p_heterogeneity, method = "BH")) %>%
    dplyr::select(metabolite, I2, p_heterogeneity, FDR)
}
