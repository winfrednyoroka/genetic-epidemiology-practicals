#!/bin/bash
#srun --nodes=1 --ntasks-per-node=14 --time=12:00:00 --mem=50gb --pty bash -i


#Updated to new version of R
module add languages/r/4.1.0

# 28.01.22 Updated to BC4 by R. Granell
# Here we use the gwas_graphs.R script to generate Q-Q and Manhattan plots
# We pass to R the UNCLEAN GWAS results file and the filename for the graphs to be saved to

workdir="$HOME/scratch/genetic-epidemiology-practicals/GWAS/scripts"
outdir="$HOME/scratch/genetic-epidemiology-practicals/GWAS/output"

Rscript ${workdir}/gwas_graphs.R ${outdir}/unclean.BMI.glm.linear.add ${outdir}/unclean_bmi
