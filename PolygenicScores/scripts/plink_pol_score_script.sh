#!/bin/bash

module load apps/plink2

cd ~/genetic-epidemiology-practicals/PolygenicScores/data

DATADIR="/newshared/genetic.epidemiology/PolScore"

plink \
	--bfile $DATADIR/data/ALSPAC_Clumped_SNPs \
	--score effect_allele_weights.txt \
	--q-score-range prs_thresholds.txt effect_allele_pvalue.txt \
	--out pgc_scz_pol_scores
     
