#!/bin/bash

# 03.01.22 Updated for BC4 by R. Granell
# Uses PLINK2

module add apps/plink/2.00
datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/data"
outdir="$HOME/scratch/genetic-epidemiology-practicals/GWAS/output"

plink \
        --bfile ${datadir}/geno_unclean \
        --glm no-x-sex \
        --pheno ${datadir}/phen.txt \
        --pheno-name BMI \
        --covar ${datadir}/covs.txt \
        --covar-name age sex \
        --out ${outdir}/unclean

#This command extracts the rows from the ADDITIVE model
awk 'NR==1 || /ADD/' ${outdir}/unclean.BMI.glm.linear | sed 's/^#CHROM/CHROM/'   > ${outdir}/unclean.BMI.glm.linear.add
