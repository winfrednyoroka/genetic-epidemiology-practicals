rm(list=ls())
setwd("~/scratch/genetic-epidemiology-practicals/Mendelian_randomization")
#setwd("~/genetic-epidemiology-practicals/Mendelian_randomization")

#These packages have already been installed for you 
#install.packages("ggplot2")
#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

# Load packages
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
# We've taken GWAS hits for BMI from the GIANT consortium, specifically the Locke Nature 2015 paper doi: 10.1038/nature14177. We defined GWAS hits as SNP-bmi associations with Pvalue<5x10^-8 (this is a conventional threshold for statistical significance in GWAS)

exp_dat <- read_exposure_data(
	filename="bmi_lock_2015.txt" ,
	sep="\t",
	phenotype_col = "trait",
	snp_col="SNP",
	beta_col="beta",
	se_col="se",
	eaf_col = "eaf",
	effect_allele_col ="effect_allele",
	other_allele_col = "other_allele",
	units_col="unit",
	samplesize_col="sample_size",
	chr_col="chr",
	pos_col="pos")

#2: Obtain the SNP-CHD summary data                  
# we've chosen to define summary genetic data for CHD using the CARDIoGRAMplusC4D Nature Genetics 2013 paper  doi: 10.1038/ng.2480.
#To create chd_deloukas_2013.txt, I first identified the SNPs associated with BMI (P<5x10^-8). I then extracted the BMI SNPs from the CARIDIoGRAMplusC4D GWAS dataset.  

out_dat <- read_outcome_data(
	filename="chd_deloukas_2013.txt" ,
	sep="\t",
	phenotype_col = "trait",
	snp_col="SNP",
	beta_col="beta",
	se_col="se",
	eaf_col = "eaf",
	effect_allele_col ="effect_allele",
	other_allele_col = "other_allele",
	units_col="unit",
	ncase_col="ncase",
	ncontrol_col="ncontrol",
	chr_col="chr",
	pos_col="pos")


#3: Harmonise the CHD and BMI datasets              

# We harmonise the CHD and BMI datasets so that the effect alleles are the same  
#Where the effect alleles are different, this syntax will flip the log odds ratio and effect alleles in the outcome dataset to reflect the effect allele in the exposure dataset

# First, we harmonise the datasets manually

dat<-merge(exp_dat,out_dat,by="SNP")
nrow(exp_dat)-nrow(dat) # notice that 8 SNPs in the exposure dataset are missing from the outcome dataset. We might lose a bit of power by dropping these SNPs 

# positions where the effect alleles are different
Pos<-which(dat$effect_allele.exposure != dat$effect_allele.outcome)
length(Pos) #we have three positions where effect alleles are different. Let's go through each one at a time

# first effect allele discrepancy
# example of non-Palindromic SNP but different reference strand 
a<-dat[Pos[1],]
a[,c("effect_allele.exposure","other_allele.exposure","eaf.exposure","effect_allele.outcome","other_allele.outcome","eaf.outcome","SNP")]
# notice that the effect alleles have been coded using different strands. If you flip the strand for the outcome dataset, is the effect alleles still different to the exposure dataset? 

strand1<-c("A","T","G","C")
strand2<-c("T","A","C","G")

# code to flip strand
ea<-strand2[strand1 %in% a$effect_allele.outcome ]
oa<-strand2[strand1 %in% a$other_allele.outcome ]
a$effect_allele.outcome<-ea
a$other_allele.outcome <-oa
a[,c("effect_allele.exposure","other_allele.exposure","eaf.exposure","effect_allele.outcome","other_allele.outcome","eaf.outcome","SNP")]
# No, if you flip the strand, then the effect alleles are the same. No further action is needed 


# second effect allele discrepancy. 
# example of Palindromic SNP with same reference strand
dat[Pos[2],c("effect_allele.exposure","other_allele.exposure","eaf.exposure","effect_allele.outcome","other_allele.outcome","eaf.outcome","SNP")]
# The effect alleles are apprently different but the SNP is also palindromic. This creates an ambiguity, in that the effect alleles could be the same but just coded using different strands. How can we resolve the ambiguity? The effect allele frequencies are very different. One is 0.32 and other is 0.67. It's clear from effect allele frequency that the effect alleles are different. We therefore need to flip the effect direction for this SNP in the outcome dataset

# Third effect allele discrepancy. 
# example of non-palindromic SNP with same reference strand

dat[Pos[3],c("effect_allele.exposure","other_allele.exposure","eaf.exposure","effect_allele.outcome","other_allele.outcome","eaf.outcome","SNP")]

# The effect alleles are clearly different between the exposure and outcome datasets.  We therefore need to flip the effect direction for this SNP in the outcome dataset 

# We now harmonise the summary data for the two discrepant SNPs where effect alleles were different between the exposure and outcome datasets

Pos<-Pos[2:3] 

ea<-dat$effect_allele.outcome[Pos]
oa<-dat$other_allele.outcome[Pos]
eaf<-dat$eaf.outcome[Pos]
beta<-dat$beta.outcome[Pos]

dat$effect_allele.outcome[Pos]<-oa
dat$other_allele.outcome[Pos]<-ea
dat$eaf.outcome[Pos]<-1-eaf
dat$beta[Pos]<-beta*-1

dat[Pos,c("effect_allele.exposure","other_allele.exposure","eaf.exposure","effect_allele.outcome","other_allele.outcome","eaf.outcome","SNP")]

# the data is harmonised for the two SNPs with discrepany effect alleles

# We can also do the harmonisation using the TwoSampleMR harmonise_data function. 

dat<-harmonise_data(exp_dat, out_dat)
Pos<-which(dat$effect_allele.exposure!= dat$effect_allele.outcome)
length(Pos) #0 rows where effect alleles are different

#4: estimate the effect of genetically elevated BMI on CHD 

# using the mr() function from the TwoSampleMR package

?mr_method_list() #useful way to see all methods available 

# Let's use MR-Base's TwoSampleMR R package to estimate the effects
res<-mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_wald_ratio","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
res # Results from the MR-Base package using various methods 

#to aid result interpretation, calculate odds ratios and 95% confidence intervals. 
res$or<-round(exp(res$b),2)
res$lci<-round(exp(res$b-1.96*res$se),2)
res$uci<-round(exp(res$b+1.96*res$se),2)
res[,c("method","or","lci","uci","pval")]

# Let's investigate what's happening "under the hood" by estimating the effects ourselves  

# multiplicative random effects IVW
# this is linear regression with intercept  constrained to pass through zero. This means that when SNP effect on exposure is zero the SNP effect on outcome is also zero. This will be true if there's no horizontal pleiotropy or if any horizontal pleioptry is "balanced"

model<-summary(lm(dat$beta.outcome ~ -1 + dat$beta.exposure,weights = 1/dat$se.outcome^2))

# see that result is same as that estimated by the mr() function in TwoSampleMR package
res[res$method == "Inverse variance weighted (multiplicative random effects)",c("b","se","pval") ]
# Notice that P values are different. 
# different P values because TwoSampleMR estimats P values using Z test while the lm function uses more conservative t test
# Z test
pnorm(abs(model$coefficients[3]),lower.tail=FALSE)*2
# t test 
pt(abs(model$coefficients[3]),df=nrow(dat)-1,lower.tail=FALSE)*2 #degrees of freedom are number of SNPs minus 1


# fixed effects IVW#


# we can estimate the fixed IVW effect using linear regression (removing the residual standard error)

#linear regression with removal of the sigma correction for heterogeneity

model$sigma #sigma is the residual standard error

se_fixed_effect_ivw<-model$coefficients[2]/model$sigma #remove the correction for the residual standard error 
# notice that this corrected standard error is same as the standard error estimated by the IVW method fixed effects method in the mr() function
se_fixed_effect_ivw 
res[res$method == "Inverse variance weighted (fixed effects)",c("b","se") ]


# We could also estimated the IVW fixed effect by meta-analysis of the Wald ratios
wr<-dat$beta.outcome/dat$beta.exposure #estimate the wald ratios
wr.se<-dat$se.outcome/dat$beta.exposure #standard errors for Wald ratios (this is equivalent to first order term of delta method and a limitation is that it doesn't take into account uncertainty in SNP exposure effect)

# IVW fixed effect:
beta_fixed_effect<-sum(wr*1/wr.se^2)/sum(1/wr.se^2) 
# standard error for IVW fixed effect:
se_fixed_effect<-sqrt(sum(1/wr.se^2)^-1)
#notice this is the same result as produced by IVW fixed effect method in mr() function:
res[res$method == "Inverse variance weighted (fixed effects)",c("b","se") ]


# notice that the fixed and multiplicative random effect IVW models estimate the same effect size but the standard error from the random effects model is larger. The standard error in the random effects model incorporates the residual standard error in the linear regression model. The residual standard error represents heterogeneity between the SNPs, which potentially reflects horizontal pleiotropy. Note that the result from random effects IVW is valid if any horizontal pleiotropy is "balanced". 

#MR Egger 
# this is linear regression without the intercept constrained to pass through zero (in contast to random and fixed effects IVW, which assume a zero intercept)

# unlike the IVW methods, we must further harmonise the exposure and outcome summary data so that the effect allele always reflects the exposure increasing allele. 

dat$beta.exposure # exposure effect sometimes negative
Pos<-which(dat$beta.exposure<0)

b.exp<-dat$beta.exposure[Pos]*-1
dat$beta.exposure[Pos]<-b.exp
b.out<-dat$beta.outcome[Pos]*-1
dat$beta.outcome[Pos]<-b.out

model<-summary(lm(dat$beta.outcome ~ dat$beta.exposure,weights = 1/dat$se.outcome^2))
model$coefficients 
# notice that this is same result generated by the mr() function for MR Egger
res[res$method == "MR Egger",c("b","se","pval") ]



#################################################################
#Visualise results and assess plausibility of MR assumptions through plots, MR-Egger intercept, tests for heterogeneity  and comparison of results across methods       #
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

#5. Create a forest plot of the results 
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw_fe","mr_ivw_mre","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
p3 <- mr_forest_plot(res_single)
p3
#save your plot using the png() function
png("forest.png")
p3
dev.off()

#6. Create a forest plot of the leave one out analysis results 
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


