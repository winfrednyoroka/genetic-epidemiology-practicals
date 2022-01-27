#!/bin/bash

# any line that starts with # is a comment line
# 27.01.22 Updated for BC4 by R. Granell

module add apps/plink/2.00
datadir="/mnt/storage/private/mrcieu/training/genetic_epidemiology/GenAssoc/data"
outdir="$HOME/scratch/genetic-epidemiology-practicals/GenAssoc/output"

	plink \
		--bfile ${datadir}/chr16 \
		--snp rs3751812 \
		--glm no-x-sex \
		--ci 0.95 \
		--pheno ${datadir}/phen.txt \
		--pheno-name ______ \
		--covar ${datadir}/covs.txt \
		--covar-name ________ \
		--out ${outdir}/_________
exit
