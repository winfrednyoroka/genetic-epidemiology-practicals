#!/bin/bash

module load apps/plink2

plink \
	--bfile ../data/geno_unclean \
	--linear \
	--pheno ../data/phen.txt \
	--pheno-name BMI \
	--covar ../data/covs.txt \
	--covar-name age sex \
	--out ../output/bmi       

awk 'NR==1 || /ADD/' ../output/bmi.assoc.linear > ../output/bmi.assoc.linear.add
