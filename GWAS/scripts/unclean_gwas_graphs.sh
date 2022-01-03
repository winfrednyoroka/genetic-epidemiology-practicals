#!/bin/bash

#export R_LIBS="~/R_libs"
#mkdir ~/R_libs

module add languages/r/4.0.3
# 03.01.22 Updated to BC4 by R. Granell
# Here we use the gwas_graphs.R script to generate Q-Q and Manhattan plots
# We pass to R the GWAS results file and the filename for the graphs to be saved to

Rscript gwas_graphs.R ~/scratch/genetic-epidemiology-practicals/GWAS/output/bmi.assoc.linear.add ~/scratch/genetic-epidemiology-practicals/GWAS/output/bmi_unclean
