rm(list=ls())
setwd("~/scratch/genetic-epidemiology-practicals/Mendelian_randomization")

# Lines 6 to 10 not needed. All these packages should be already installed.
# Packages needed for this practical 
install.packages("plyr")
install.packages("ggplot2")

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

# Load packages
library(plyr) #for the rbind.fill function
library(ggplot2)
library(TwoSampleMR)

#################################################################################################
#Example: MR study of the effect of body mass index on coronary heart  disease                                              
#################################################################################################
# Coronary heart disease is one of the leading causes of death and morbidity in the world. Of the known risk factors, body mass index (BMI) is one of the most important. In prospective observational studies (Wormser et al Lancet 2011 377:1085), each standard deviation (SD) (~4.56 kg/m2) increase in BMI is associated with a relative risk of 1.23 (95% CI: 1.17–1.29) for coronary heart disease (CHD), after adjustment for conventional vascular risk factors, including smoking, cholesterol, diabetes, and blood pressure. However, whether the association reflects a causal effect of BMI on CHD is unknown.
#########################################################
#Estimate the effect of genetically elevated BMI on CHD  
#########################################################

# 1: Define a genetic instrument for BMI 
ao<-available_outcomes(access_token=NULL)
# What does the available data look like? 
head(ao)

# Select the BMI data from MR-Base 
 #There are many datasets to choose from. We will choose the BMI dataset generated in the largest analysis 
ao.bmi<-ao[ao$trait=="Body mass index",]
ao.bmi<-ao.bmi[order(ao.bmi$sample_size,decreasing=TRUE),]
ao.bmi<-ao.bmi[!duplicated(ao.bmi$trait),]
data.frame(ao.bmi) #the chosen dataset correspodns to Lock et al 2015 

# note the unit of the SNP effet is SD units
ao.bmi$unit

# get the study ID for BMI
id.bmi<-ao.bmi$id

# Extract BMI instruments 
bmi_exp_data <- extract_instruments(outcomes=id.bmi,access_token=NULL)

#2: Obtain the SNP-CHD summary data                  

#There are many datasets to choose from. We will choose the CHD dataset with the largest number of cases  
ao.chd<-ao[ao$trait=="Coronary heart disease",]
ao.chd<-ao.chd[order(ao.chd$ncase,decreasing=TRUE),]
ao.chd<-ao.chd[!duplicated(ao.chd$trait),]
data.frame(ao.chd) #we've chosen the CARDIoGRAMplusC4D analysis from 2013
id.chd<-ao.chd$id
chd_out_data<-extract_outcome_data(bmi_exp_data$SNP,id.chd,proxies=1,access_token=NULL)    

#3: Harmonise the CHD and BMI datasets              

# We harmonise the CHD and BMI datasets so that the effect alleles are the same  
# This syntax will flip the log odds ratio and effect alleles in the CARDIoGRAM dataset where the effect alleles are different between CARDIoGRAMplusC4D and GIANT
dat <- harmonise_data(bmi_exp_data, chd_out_data)

# If you explore the dataset you'll notice that effect alleles and log odds ratios have been flipped in the CHD dataset where the effect allele in the CHD dataset was different from the effect allele in the BMI dataset
head(dat)



#4: estimate the effect of genetically elevated BMI on CHD

?mr_method_list() #useful way to see all methods available 

# Let's use MR-Base's TwoSampleMR R package to estimate the effects
res<-mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_wald_ratio","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
res # Results from the MR-Base package using various methods 

# calculate  odds ratios and 95% confidence intervals for the odds ratios
res$or<-round(exp(res$b),2)
res$lci<-round(exp(res$b-1.96*res$se),2)
res$uci<-round(exp(res$b+1.96*res$se),2)

res[,c("method","or","lci","uci")]

#5. Let's estimate the Wald ratio for each SNP 
res_single <- mr_singlesnp(dat)
res_single

#################################################################
#Visualise results and assess plausibility of MR assumptions through plots, MR-Egger intercept and tests for heterogeneity         #
#################################################################

#1. MR-Egger intercept test 

egg.int<-mr_pleiotropy_test(dat) 

egg_int<-round(egg.int$egger_intercept,3)
egg_lci<-round(egg.int$egger_intercept-1.96*egg.int$se,3)
egg_uci<-round(egg.int$egger_intercept+1.96*egg.int$se,3)

#When the SNP effect on BMI is zero (the intercept), then the SNP effect on CHD should also be zero. If it isn't, this is evidence for the presence of horizontal pleiotropy. The intercept can also be thought of as an estimate of the average (horizontally) pleiotropic effect. If the instruments are valid, or horizontal pleiotrpy is balanced, it should be 0. We also have to assume the InSIDE assumption (that the horizontal pleiotropic effects are not correlated with the SNP-exposure effects)  
c(egg_int,egg_lci,egg_uci)

#2. Cochran's Q test for heterogeneity 
mr_heterogeneity(dat)

# 3. Create a scatter plot of the SNP-CHD and SNP-BMI associations. Does the SNP-CHD association increase linearly as the SNP-BMI association increases? What could deviations from linearity mean? Are there any unsual data points?
p1 <- mr_scatter_plot(res, dat)
p1 #have a look at the scatter plot
#save your plot using the png() function
png("scatter.png")
p1
dev.off()

# 4.	Create a funnel plot of the results. Does the funnel plot look symmetric? What could assymetry mean? Are there any outliers?
res_single <- mr_singlesnp(dat)
p2 <- mr_funnel_plot(res_single)
p2
#save your plot using the png() function
png("funnel.png")
p2
dev.off()

# 5. Create a forest plot of the results 
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw_fe","mr_ivw_mre","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
p3 <- mr_forest_plot(res_single)
p3
#save your plot using the png() function
png("forest.png")
p3
dev.off()

# 4. Create a forest plot of the leave one out analysis results 
res_loo <- mr_leaveoneout(dat)
p4 <- mr_leaveoneout_plot(res_loo)
p4
png("loo.png")
p4
dev.off()

#################################################################
#results interpretation 
#################################################################
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


