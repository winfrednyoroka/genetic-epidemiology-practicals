#!/bin/bash

# 28.01.22 Updated for BC4 by R. Granell
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
        --out ${outdir}/clean

# Extract rows from ADDITIVE model and remove # from header
awk 'NR==1 || /ADD/' ${outdir}/clean.BMI.glm.linear | sed 's/^#CHROM/CHROM/' > ${outdir}/clean.BMI.glm.linear.add

# Extract chromosome 16 data for Exercise 4 (locus plot)
awk 'NR==1 || $1==16' ${outdir}/clean.BMI.glm.linear.add > ${outdir}/clean16.BMI.glm.linear.add
