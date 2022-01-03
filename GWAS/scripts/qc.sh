#!/bin/bash

# 03.01.22 Updated for BC4 by R. Granell
# Uses PLINK2

module add apps/plink/2.00
datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GWAS/data"

plink \
 	--bfile ${datadir}/geno_unclean \
 	--maf 0.01 \
 	--hwe 1e-6 \
 	--geno 0.05 \
 	--mind 0.05 \
 	--make-bed \
 	--out ~/scratch/genetic-epidemiology-practicals/GWAS/data/geno_qc

