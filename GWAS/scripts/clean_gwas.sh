#!/bin/bash

# 03.01.22 Updated for BC4 by R. Granell
# USes PLIN2

module add apps/plink/2.00
datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/data"

plink \
        --bfile ~/scratch/genetic-epidemiology-practicals/GWAS/data/geno_qc \
        --glm no-x-sex \
        --pheno ${datadir}/phen.txt \
        --pheno-name BMI \
        --covar ${datadir}/covs.txt \
        --out ~/scratch/genetic-epidemiology-practicals/GWAS/output/bmi_clean

awk 'NR==1 || /ADD/' ~/scratch/genetic-epidemiology-practicals/GWAS/output/bmi_clean.BMI.glm.linear > ~/scratch/genetic-epidemiology-practicals/GWAS/output/bmi_clean.BMI.glm.linear.add
