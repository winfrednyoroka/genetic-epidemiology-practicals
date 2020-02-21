#!/bin/bash

module load apps/plink2
datadir="/newshared/genetic.epidemiology/GWAS/data"

plink \
	--bfile ${datadir}/geno_unclean \
	--linear \
	--pheno ${datadir}/phen.txt \
	--pheno-name BMI \
	--covar ${datadir}/covs.txt \
	--covar-name age sex \
	--out ~/genetic-epidemiology-practicals/GWAS/output/bmi       

awk 'NR==1 || /ADD/' ~/genetic-epidemiology-practicals/GWAS/output/bmi.assoc.linear > ~/genetic-epidemiology-practicals/GWAS/output/bmi.assoc.linear.add
