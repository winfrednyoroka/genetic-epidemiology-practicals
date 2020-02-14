#!/bin/bash

# any line that starts with # is a comment line
module load apps/plink2

	plink \
		--bfile ~/pract2_PopgenAssoc/data/chr16 \
		--snp rs3751812 \
		--linear \
		--ci 0.95 \
		--pheno ~/pract2_PopgenAssoc/data/phen.txt \
		--pheno-name BMI \
		--covar ~/pract2_PopgenAssoc/data/covs.txt \
		--covar-name sex,age \
		--out ~/pract2_PopgenAssoc/output/rs3751812_BMI
