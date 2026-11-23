## === R/01_plsda_paired_ttest.R ===
## Stage-specific metabolite identification via PLS-DA and paired t-test
##
## Input:  serum_metabolites — data.frame [id, trimester, metabolite_1 ... metabolite_n]
## Output: log2FC, VIP, P value, FDR for each metabolite per transition

library(dplyr)
library(ropls)

perform_analysis <- function(data, time1, time2, seed) {
  data_raw <- data %>% filter(trimester %in% c(time1, time2))
  data_log <- data_raw %>% mutate(across(-c(id, trimester), ~ scale(log2(.)) %>% as.vector()))
  
  raw_t1 <- data_raw %>% filter(trimester == time1) %>% arrange(id)
  raw_t2 <- data_raw %>% filter(trimester == time2) %>% arrange(id)
  log_t1 <- data_log %>% filter(trimester == time1) %>% arrange(id)
  log_t2 <- data_log %>% filter(trimester == time2) %>% arrange(id)
  metabolites <- setdiff(colnames(raw_t1), c("id", "trimester"))
  
  fc   <- sapply(metabolites, function(x) mean(raw_t2[[x]] / raw_t1[[x]], na.rm = TRUE))
  pval <- sapply(metabolites, function(x) t.test(log_t1[[x]], log_t2[[x]], paired = TRUE)$p.value)
  fdr  <- p.adjust(pval, method = "BH")
  
  X <- as.matrix(rbind(log_t1[, metabolites], log_t2[, metabolites]))
  Y <- factor(c(rep(time1, nrow(log_t1)), rep(time2, nrow(log_t2))))
  set.seed(seed)
  vip <- getVipVn(opls(X, Y))
  
  data.frame(metabolite = metabolites, log2FC = log2(fc), VIP = vip,
             p_value = pval, FDR = fdr, transition = paste0(time1, "_to_", time2))
}

# Individual log2FC adjusted for sampling interval and baseline intensity
metabolite_change <- function(earlier, later, interval) {
  log2fc <- log2(later / earlier)
  residuals(lm(log2fc ~ interval + earlier, na.action = na.exclude))
}

# Overall distance between matched trimester profiles
compare_distances <- function(t1, t2, t3) {
  z1 <- scale(log2(as.matrix(t1)))
  z2 <- scale(log2(as.matrix(t2)))
  z3 <- scale(log2(as.matrix(t3)))
  d1 <- sqrt(rowSums((z2 - z1)^2))
  d2 <- sqrt(rowSums((z3 - z2)^2))
  list(distances = data.frame(T1_T2 = d1, T2_T3 = d2),
       test = wilcox.test(d1, d2, paired = TRUE))
}
