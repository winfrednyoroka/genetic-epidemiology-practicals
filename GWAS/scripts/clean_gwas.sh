#!/bin/bash

# 03.01.22 Updated for BC4 by R. Granell
# USes PLIN2

module add apps/plink/2.00
datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/data"
outdir="$HOME/scratch/genetic-epidemiology-practicals/GWAS/output"

plink \
        --bfile ~/scratch/genetic-epidemiology-practicals/GWAS/data/geno_qc \
        --glm no-x-sex \
        --pheno ${datadir}/phen.txt \
        --pheno-name BMI \
        --covar ${datadir}/covs.txt \
        --out ${outdir}/bmi_clean

# Extract rows from ADDITIVE model and remove # from header
awk 'NR==1 || /ADD/' ${outdir}/bmi_clean.BMI.glm.linear | sed 's/^#CHROM/CHROM/' > ${outdir}/bmi_clean.BMI.glm.linear.add
awk 'NR==1 || $1==16' ${outdir}/bmi_clean.BMI.glm.linear.add > ${outdir}/bmi_clean16.BMI.glm.linear.add
