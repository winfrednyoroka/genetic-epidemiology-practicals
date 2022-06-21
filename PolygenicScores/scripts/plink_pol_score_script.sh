#!/bin/bash

module load apps/plink2

cd ~/scratch/genetic-epidemiology-practicals/PolygenicScores/data

DATADIR="/mnt/storage/private/mrcieu/training/genetic_epidemiology/PolScore"

plink \
	--bfile $DATADIR/data/ALSPAC_Clumped_SNPs \
	--score effect_allele_weights.txt \
	--q-score-range prs_thresholds.txt effect_allele_pvalue.txt \
	--out pgc_scz_pol_scores
