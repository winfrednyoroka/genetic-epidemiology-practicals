################################################################################################################
# Step 1: Generate files needed to construct polygenic scores using discovery sample summary statistics
################################################################################################################

# Change your working directory (you may need to change the filepath in next line)
setwd("~/scratch/genetic-epidemiology-practicals/PolygenicScores/data")

# Open discovery sample file with 180,290 SNPs from PGC schizophrenia meta-analysis
SCZdata <- read.table("pgc_scz_gwas_summary_statistics.txt", header=T)
# Explore the PGC summary stats
names(SCZdata)
head(SCZdata)
dim(SCZdata)

# As schizophrenia was modelled as a binary outcome (i.e. cases versus controls) in the PGC GWAS, the estimated 
# effects for each SNP are odds ratios (ORs). To generate weights for each effect allele (A1) we must first 
# convert the ORs into log(ORs):
SCZdata$weight <- log(SCZdata$OR)


# Next, you will need to generate 3 files that are needed for constructing the schizophrenia polygenic scores in PLINK 1.90

# A) effect_allele_weights.txt (input for plink --score command). This is a tab-delimited file containing the SNP IDs, 
# effects alleles and corresponding weights. To generate this file type:
write.table(subset(SCZdata, select=c("SNP","A1","weight")), file="effect_allele_weights.txt", 
            row.names=F, col.names=F, quote=F, sep ="\t")

#	B) effect_allele_pvalue.txt (input for plink --q-score-range command). This is a tab-delimited file containing the SNP ID 
# and their corresponding p-values (for the association between the SNPs and schizophrenia). To generate this file type:
write.table(subset(SCZdata, select=c("SNP","P")), file="effect_allele_pvalue.txt", 
            row.names=F, col.names=F, quote=F, sep ="\t")  


################################################################################################################
# Step 3: Prepare phenotype file in target sample 
################################################################################################################

# Change your working directory (you may need to change the filepath in next line)
setwd("~/scratch/genetic-epidemiology-practicals/PolygenicScores/data")

# Explore the phenotypic data (this data set includes a BMI variable as well as blood pressure, 
# hypertension and C-reactive protein measures):
outcomeData <- read.table("outcome.txt", header=T)

# Examine the distribution of the BMI variable
summary(outcomeData$bmi)
# visualise the distribution as a histogram and save as jpeg
jpeg(file="~/scratch/genetic-epidemiology-practicals/PolygenicScores/results/hist_bmi.jpeg")
hist(outcomeData$bmi, breaks=40, xlim=c(0,max(outcomeData$bmi)), col="lightgrey")
dev.off()

# Replace the BMI outliers (BMI greater than 100 or less than 10) with missing values (NA)
outcomeData$bmi[outcomeData$bmi>100 | outcomeData$bmi<10] <- NA

# Examine the distribution of the BMI variable after removing the outliers: 
summary(outcomeData$bmi)
# Again, visualise the distribution as a histogram and save as jpeg
jpeg(file="~/scratch/genetic-epidemiology-practicals/PolygenicScores/results/hist_bmi_clean.jpeg")
hist(outcomeData$bmi, breaks=40, xlim=c(0,max(outcomeData$bmi, na.rm=TRUE)), col="lightgrey")
dev.off()

# save the cleaned phenotypic data
write.table(outcomeData, file="outcome_clean.txt", row.names=F, col.names=T, quote=F, sep ="\t")  




