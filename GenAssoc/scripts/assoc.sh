#!/bin/bash

# any line that starts with # is a comment line
module load apps/plink2
datadir="/newshared/genetic.epidemiology/GenAssoc/data"

	plink \
		--bfile ${datadir}/chr16 \
		--snp rs3751812 \
		--linear \
		--ci 0.95 \
		--pheno ${datadir}/phen.txt \
		--pheno-name BMI \
		--covar ${datadir}/covs.txt \
		--covar-name sex,age \
		--out ~/genetic-epidemiology-practicals/GenAssoc/output/rs3751812_BMI
