rm(list=ls())
setwd("~/MR-practice/")

# Packages needed for this practical 
install.packages("devtools")
install.packages("plyr")
install.packages("ggplot2")
install_github("MRCIEU/TwoSampleMR@0.4.25")

# Load packages
library(devtools)
library(plyr) #for the rbind.fill function
library(ggplot2)
library(TwoSampleMR)

#################################################################################################
#Practice 1: One exposure on one outcome                                                        #              
#################################################################################################

#########################################################
#Stage 1: Select SNP-exposure summary data              #              
#########################################################
rm(list=ls())
# You can extract summary data for several traits directly from MR-base (for details, check: https://mrcieu.github.io/TwoSampleMR/)
ao<-available_outcomes(access_token=NULL)
# What does the available data look like? 
head(ao)

# Select the BMI data from MR-Base database 
ao.bmi<-ao[ao$trait=="Body mass index",]
ao.bmi<-ao.bmi[ao.bmi$year==2015,]
ao.bmi<-ao.bmi[ao.bmi$population=="Mixed",]
id.bmi<-ao.bmi$id

# Extract BMI instruments 
bmi_exp_data <- extract_instruments(outcomes=id.bmi,access_token=NULL)

##################################################################################################
# Stage 2: Select SNP-outcome summary data and harmonise datasets                                #
##################################################################################################

# Select the CHD data from MR-Base database 
ao.chd<-ao[ao$trait=="Coronary heart disease",]
ao.chd<-ao.chd[ao.chd$year==2015,]
id.chd<-ao.chd$id
chd_out_data<-extract_outcome_data(bmi_exp_data$SNP,id.chd,proxies=1,access_token=NULL)    

# d. Harmonise the CHD and BMI datasets so that the effect alleles are the same (and reflect the BMI increasing allele) 
# This syntax will flip the log odds ratio and effect alleles in the CARDIoGRAM dataset where the effect alleles are different between CARDIoGRAMplusC4D and GIANT
dat <- harmonise_data(bmi_exp_data, chd_out_data)

# If you explore the dataset you'll notice that effect alleles and log odds ratios have been flipped in the CHD dataset where the effect allele in the CHD dataset was different from the effect allele in the BMI dataset
head(dat)


#################################################################
#Stage 3: estimate the effect of genetically elevated BMI on CHD#
#################################################################

# Let's use the MR-Base R package to estimate the effects
res<-mr(dat,method_list=c("mr_ivw","mr_wald_ratio","mr_egger_regression","mr_weighted_median"))
res # Results from the MR-Base package using various methods including MR-Egger and weighted median sensitivity analyses   
#estimate odds ratio and 95% confidence interval
exp(res$b[1])
exp(res$b[1]-1.96*res$se[1])
exp(res$b[1]+1.96*res$se[1])

egg.int<-mr_pleiotropy_test(dat) # MR-Egger intercept test 
egg.int

# Let's estimate the Wald ratio for each SNP 
res_single <- mr_singlesnp(dat)
res_single

#################################################################
#Stage 4: visulise the causal estimate for BMI on CHD           #
#################################################################

# 1. Create a scatter plot of the SNP-CHD and SNP-BMI associations. Does the SNP-CHD association increase linearly as the SNP-BMI association increases? What could deviations from linearity mean? Are there any unsual data points?
p1 <- mr_scatter_plot(res, dat)
p1 #have a look at the scatter plot
#save your plot using the png() function
png("scatter.png")
p1
dev.off()

# 2.	Create a funnel plot of the results. Does the funnel plot look symmetric? What could assymetry mean? Are there any outliers?
res_single <- mr_singlesnp(dat)
p2 <- mr_funnel_plot(res_single)
p2
#save your plot using the png() function
png("funnel.png")
p2
dev.off()

# 3. Create a forest plot of the results 
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_egger_regression","mr_weighted_median"))
p3 <- mr_forest_plot(res_single)
p3
#save your plot using the png() function
png("forest.png")
p3
dev.off()

# 3. Create a forest plot of the leave one out analysis results 
res_loo <- mr_leaveoneout(dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4
png("loo.png")
p4
dev.off()


# 	Do you think BMI causes CHD? Split into 3 groups and discuss. One individual from each group will present a summary of the discussion to the class. Consider these questions in your discussion:
# a.	What is the odds ratio for coronary heart disease per unit increase in genetically elevated BMI? 
# b.	Is there evidence for pleiotropy or violations of MR assumptions? 
# c.	Are the genetic and observational effects directionally similar? Are they comparable? Note: In prospective observational studies, the multivariate adjusted relative risk for coronary heart disease per standard deviation higher BMI is 1.23 (95% CI: 1.17–1.29) (Wormser et al Lancet 2011 377:1085). 
# d.	Are these results compatible with a causal effect of BMI on coronary heart disease?
# e.	Can you think of reasons for caution?

# Also consider these assumptions and limitations of Mendelian randomization in your interpretation of the results: 

# Core assumptions of Mendelian randomization: 
# 1) The instrument is associated with the exposure
# 2) The instrument is not associated with confounders of the exposure-outcome association
# 3) The instrument is associated with the outcome exclusively via its effect on the exposure (also known as the exclusion restriction assumption)

# Other assumptions and limitations
# 1)	Samples should be independent 
# 2)	Two samples should be from the same population 
# 3)	Potential impact of weak instruments bias 
# 4)	Potential impact of winner’s curse, whereby the effects of the SNPs on BMI may be overestimated in the discovery stage of the GWAS 
# 5)	Strand ambiguity (when two GWAS studies use different reference strands)
# 6)	SNPs should be independent 


