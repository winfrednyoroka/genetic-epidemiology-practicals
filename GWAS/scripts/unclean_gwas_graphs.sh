#!/bin/bash

export R_LIBS="~/R_libs"
mkdir ~/R_libs
module load languages/R-3.6.1-gcc9.1.0

# Here we use the gwas_graphs.R script to generate Q-Q and Manhattan plots
# We pass to R the GWAS results file and the filename for the graphs to be saved to

Rscript gwas_graphs.R ~/genetic-epidemiology-practicals/GWAS/output/bmi.assoc.linear.add ~/genetic-epidemiology-practicals/GWAS/output/bmi_unclean
