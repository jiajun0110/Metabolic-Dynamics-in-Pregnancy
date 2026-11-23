#!/bin/bash
## GWAS of metabolite changes and LD clumping
## Covariates: age, pre-pregnancy BMI, parity, gravidity and PC1-PC10.

set -e
BFILE="data/genotype"
PHENO_DIR="data/gwas_phenotypes"
OUT_DIR="output/gwas"
mkdir -p "$OUT_DIR"

# GCTA fastGWA-lr
for i in $(seq 1 38); do
  gcta64 --bfile "$BFILE" --fastGWA-lr --autosome \
    --pheno "${PHENO_DIR}/pheno_${i}.txt" \
    --qcovar "${PHENO_DIR}/qcovar.txt" \
    --out "${OUT_DIR}/dynamics_${i}"
done

# LD clumping
RESULT="${OUT_DIR}/clumped_results.txt"
printf 'Trait\tCHR\tPOS\tSNP\tBETA\tSE\tP\tAF1\tN\n' > "$RESULT"
for i in $(seq 1 38); do
  GWAS="${OUT_DIR}/dynamics_${i}.fastGWA"
  rm -f "${OUT_DIR}/clump_${i}.clumped"
  plink --bfile "$BFILE" --clump "$GWAS" \
    --clump-p1 5e-8 --clump-p2 1 --clump-r2 0.05 --clump-kb 500 \
    --clump-snp-field SNP --clump-field P --out "${OUT_DIR}/clump_${i}"
  if [ -f "${OUT_DIR}/clump_${i}.clumped" ]; then
    awk -v trait="$i" 'BEGIN {OFS="\t"}
      NR==FNR {if(FNR>1) row[$2]=trait OFS $1 OFS $3 OFS $2 OFS $8 OFS $9 OFS $10 OFS $7 OFS $6; next}
      FNR>1 && NF>0 {if($3 in row) print row[$3]}
    ' "$GWAS" "${OUT_DIR}/clump_${i}.clumped" >> "$RESULT"
  fi
done

# Study-wide significance
awk 'NR==1 || $7 < 5e-8/38' "$RESULT" > "${OUT_DIR}/study_wide_results.txt"
