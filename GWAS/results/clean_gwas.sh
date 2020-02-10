#!/bin/bash

module load apps/plink2

plink \
	--bfile ../output/geno_qc \
	--linear \
	--pheno ../data/phen.txt \
	--pheno-name BMI \
	--covar ../data/covs.txt \
	--out ../output/bmi_clean   


awk 'NR==1 || /ADD/' ../results/bmi_clean.assoc.linear > ../output/bmi_clean.assoc.linear.add
