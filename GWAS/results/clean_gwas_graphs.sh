#!/bin/bash

export R_LIBS="~/R_libs"
mkdir ~/R_libs
module load languages/R-3.2.2-ATLAS

# Here we use the gwas_graphs.R script to generate Q-Q and Manhattan plots
# We pass to R the GWAS results file and the filename for the graphs to be saved to

Rscript ../scripts/gwas_graphs.R ../results/bmi_clean.assoc.linear.add ../output/bmi_clean
