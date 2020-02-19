#!/bin/bash

module load apps/plink2

plink \
 	--bfile ${datadir}/geno_unclean \
 	--maf 0.01 \
 	--hwe 1e-6 \
 	--geno 0.05 \
 	--mind 0.05 \
 	--make-bed \
 	--out ~/genetic-epidemiology-practicals/GWAS/data/geno_qc

