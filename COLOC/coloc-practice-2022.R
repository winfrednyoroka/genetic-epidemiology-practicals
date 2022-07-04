###########################################
#practice for COLOC analysis - Chris Zheng#
###########################################

#########Approximate Bayes Factor colocalisation analyses############
#H0: neither trait has a genetic association in the region          #
#H1:only trait 1 has a genetic association in the region            #
#H2: only trait 2 has a genetic associaiton in the region           #
#H3: both traits are associated, but with different causal variants #
#H4: both traits are associated, and share a single causal variant  #
#####################################################################

################################Parameters for coloc#####################################
#beta1 - effect sizes for trait 1                                                       #
#beta2 - effect sizes for trait 2                                                       #
#p1 - p-values for trait 1                                                              #
#p2 - p-values for trait 2                                                              #
#MAF1 - minor allele frequencies for SNPs of trait 1                                    #
#MAF2 - minor allele frequencies for SNPs of trait 2                                    #
#N1 - sample size for trait 1                                                           #
#N2 - sample size for trait 2                                                           #
#s - for a case control dataset, the proportion of samples in dataset that are cases  #
#########################################################################################


##install and call packages (these packages are already been installed for you)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("snpStats")


#install.packages("coloc")
#devtools::install_github("chr1swallace/coloc") #alternative approach

library(coloc)
#library(ggplot2)
rm(list=ls(all=TRUE)) 

#change working directory

#Set working directory to where the coloc input files are stored.
#setwd("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/Teaching/Genetic-epidemiology/2018/practice/COLOC")
setwd("~/scratch/genetic-epidemiology-practicals/COLOC")
################## ANALYSIS 1: colocalization analysis of exposure level of CSF1 on years of schooling ##################
##coloc function for quantitative trait## 
coloc.analysis.quant <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2){
  
  #Convert the inputs in order to run in coloc function.
  #type, quant (quantitative) for pQTL study 
  dataset1 <- list(beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  
  
  #type, quant (quantitative) for quantitative trait study 
  dataset2 <- list(beta=beta2, varbeta=se2^2,MAF=MAF2, type="quant",N=N2)
  
  
  #Run the coloc analysis, setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12) as 1E-5.
  #p1 prior probability a SNP is associated with trait 1, default 1e-4
  #p2 prior probability a SNP is associated with trait 2, default 1e-4
  #p12 prior probability a SNP is associated with both traits, default 1e-5
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  
  #Format the data to save out.
  
  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  
  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  
  #Make the filename and save out.
  return(df)
  
}

#call the input file
input_file <- "CSF1-Schooling.txt"
df <- read.table(input_file, header=T)
df <- df[complete.cases(df), ]

#Get the MAF.  
df$MAF1 <- NULL
df$MAF2 <- NULL
for (j in 1:nrow(df)){
    if(df$eaf.exposure[j]>0.5){df$MAF1[j]=1-df$eaf.exposure[j]} else {df$MAF1[j]=df$eaf.exposure[j]}  
    if(df$eaf.outcome[j]>0.5){df$MAF2[j]=1-df$eaf.outcome[j]} else {df$MAF2[j]=df$eaf.outcome[j]}  
}
    
#sample size for two studies
N1<-3300
N2<-293723
      
result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, N1, N2) 
result


################## ANALYSIS 2: colocalization analysis of protein level of MUC16 and asthma ##################
##coloc function for binary trait##
coloc.analysis <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,s){
  
  #Convert the inputs in order to run in coloc function.
  #type, quant (quantitative) for pQTL study
  dataset1 <- list(beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  
  #type, cc (case coontrol) for binary study
  dataset2 <- list(beta=beta2, varbeta=se2^2,MAF=MAF2, type="cc",s=s,N=N2)
  
  #Run the coloc analysis, setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12) as 1E-5.
  #p1 prior probability a SNP is associated with trait 1, default 1e-4
  #p2 prior probability a SNP is associated with trait 2, default 1e-4
  #p12 prior probability a SNP is associated with both traits, default 1e-5
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  
  #Format the data to save out.
  
  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  
  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  
  #Make the filename and save out.
  return(df)

}

#call the input file
input_file <- "MUC16-Asthma.txt"
df <- read.table(input_file, header=T)
df <- df[complete.cases(df), ]

#Get the MAF.
df$MAF1 <- NULL
df$MAF2 <- NULL
for (j in 1:nrow(df)){
  if(df$eaf.exposure[j]>0.5){df$MAF1[j]=1-df$eaf.exposure[j]} else {df$MAF1[j]=df$eaf.exposure[j]}  
  if(df$eaf.outcome[j]>0.5){df$MAF2[j]=1-df$eaf.outcome[j]} else {df$MAF2[j]=df$eaf.outcome[j]}  
}

#sample size for two studies 
N1<-3300
N2<-293723
#the proportion of samples in study 2 that are cases. Number of cases = 39049; total sample size = 293723.
s<- 39049/293723
  
result2 <- coloc.analysis(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, N1, N2, s) 
result2

