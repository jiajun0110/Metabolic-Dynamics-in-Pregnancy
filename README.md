# Metabolic dynamics across pregnancy and associations with maternal health and fetal growth

Analysis code for: **Zhao et al. Metabolic dynamics across pregnancy and associations with maternal health and fetal growth.**

## Overview

Code for characterizing stage-specific metabolic dynamics across gestation in 1,321 pregnancies from four cohorts (THSBC, HBC, WeBirth, DPC), and their associations with maternal and fetal health outcomes and postnatal weight.

## Repository structure

```
R/
├── 01_plsda_paired_ttest.R               # Stage-specific metabolite identification
├── 02_I2_heterogeneity.R                 # I² statistics for stage-exclusivity
├── 03_network_analysis.R                 # Intra-stage co-variation networks
├── 04_cross_stage_permanova.R            # Cross-stage dependency (PERMANOVA)
├── 05_kmeans_GMT.R                       # K-means clustering for GMT types
├── 06_lasso_prediction.R                 # LASSO + XGBoost outcome prediction
├── 07_binary_classification.R            # Logistic models for GDM/PIH/LBW/macrosomia
├── 08_fetal_growth_trajectory.R          # Linear mixed models for fetal growth
├── 09_postnatal_weight.R                # Postnatal weight models
├── 10_permanova_variance_decomposition.R # Backward-selection PERMANOVA (multi-omics)
├── 11_lasso_per_metabolite_decomposition.R # Per-metabolite adjusted R² decomposition
├── 13_species_metacyc_associations.R     # Species and MetaCyc associations
├── 14_bmi_mediation.R                    # BMI-dynamics mediation analysis
bash/
├── 12_gwas.sh                            # GCTA fastGWA + LD clumping
```

## Scripts and corresponding analyses

| Script | Analysis |
|--------|----------|
| 01 | PLS-DA (VIP) + paired t-test for stage-specific metabolites |
| 02 | I² heterogeneity confirming temporal specificity |
| 03 | Spearman correlation networks, degree centrality, clustering coefficients |
| 04 | Per-metabolite PERMANOVA for cross-stage dependency |
| 05 | K-means (NbClust) + Jaccard bootstrap for GMT-1/GMT-2 |
| 06 | LASSO + XGBoost, SHAP and sampling-time sensitivity | 
| 07 | Logistic models, disease-specific CRFs, bootstrap AUC comparison and GMT associations |
| 08 | LMM for fetal weekly growth rate; partial correlations; tertiles and FGR |
| 09 | WAZ LMM, weight trajectories and monthly weight gain |
| 10 | Backward-selection PERMANOVA across genetics, microbiome, maternal factors, diet |
| 11 | LASSO for non-genetic factors, GWAS-selected variants and adjusted R² decomposition |
| 12 | GCTA fastGWA-lr for 38 dynamics phenotypes + LD clumping |
| 13 | Species–metabolite and MetaCyc associations |
| 14 | BMI–metabolite association and mediation of the BMI–OGTT2h association |

## Software

### R (v4.4.1)

| Package | Version | Purpose |
|---------|---------|---------|
| ropls | 1.38.0 | PLS-DA and VIP |
| glmnet | 4.1.8 | LASSO |
| xgboost | 1.7.11.1 | Outcome prediction and SHAP |
| lme4 / lmerTest | 1.1.37 | Linear mixed models |
| vegan | 2.6.10 | PERMANOVA |
| igraph | 2.1.4 | Network metrics |
| Hmisc | 5.2.3 | Correlation matrices |
| NbClust | 3.0.1 | Optimal cluster number |
| fpc | 2.2.13 | Jaccard bootstrap stability |
| mediation | 4.5.0 | Causal mediation |
| pROC | 1.18.5 | AUC |

Additional: dplyr, tidyr, stringr, purrr, broom, boot, data.table, tibble, ppcor

### Command-line tools

| Tool | Version | Purpose |
|------|---------|---------|
| GCTA | 1.94.1 | fastGWA-lr |
| PLINK | 1.90b6.7 | LD clumping |
| PLINK2 | 2.0.0-a.6.4LM | GWAS for MR instruments |

## Notes

- Scripts present statistical analysis logic only. Data loading, column renaming, and visualization code are omitted.
- Input data formats are described in the header comments of each script.
- Users should adapt variable names and file paths to their own data structures.

## Data availability

Metagenomic sequencing and genotype data have been deposited as described in the manuscript. Processed data required to reproduce the analyses are available upon reasonable request.

## Contact

zhaojiajun@westlake.edu.cn
