#!/bin/bash
#run --nodes=1 --ntasks-per-node=14 --time=12:00:00 --mem=50gb --pty bash -i

#export R_LIBS="~/R_libs"
#mkdir ~/R_libs

module add languages/r/4.0.3
#module add X11/20171023-GCCcore-6.4.0
#module load X11/20171023-GCCcore-6.4.0

# 03.01.22 Updated to BC4 by R. Granell
# Here we use the gwas_graphs.R script to generate Q-Q and Manhattan plots
# We pass to R the GWAS results file and the filename for the graphs to be saved to

outdir="$HOME/scratch/genetic-epidemiology-practicals/GWAS/output"

xvfb-run --server-args="-screen 0 1024x768x24" Rscript gwas_graphs.R ${outdir}/bmi.BMI.glm.linear.add ${outdir}/bmi_unclean
#Rscript gwas_graphs.R ${outdir}/bmi.BMI.glm.linear.add ${outdir}/bmi_unclean
