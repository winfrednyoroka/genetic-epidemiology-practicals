#!/bin/bash

# any line that starts with # is a comment line
module load apps/plink2

	plink \
		--bfile ${datadir}/chr16 \
		--snp rs3751812 \
		--linear \
		--ci 0.95 \
		--pheno ${datadir}/phen.txt \
		--pheno-name ______ \
		--covar ${datadir}/covs.txt \
		--covar-name ______ \
		--out ~/genetic-epidemiology-practicals/GenAssoc/output/rs3751812_BMI
