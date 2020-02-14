# Change your working directory (you may need to change the filepath in next line)
setwd("E:/GenEpiSC_2020/17_PRS_PRACT/data")

# Remove everything in R workspace
rm(list=ls())

# Load in cleaned outcome data
outcomeData <- read.table("outcome_clean.txt", header=T)

# Load in each PRS file, rename columns, generate standardised score and drop PHENO and FID columns
for(i in 1:12) {
  infile <- paste0("pgc_scz_pol_scores.S", i,".profile",sep="")
  assign(paste0("S",i), read.table(infile, header=T))
  d<-get(paste0("S",i))
  scorenam <- paste0("SCORE_S", i)
  scorestd <- paste0("stdSCORE_S", i)
  cntnam <- paste0("CNT_S", i)
  cnt2nam <- paste0("CNT2_S", i)
  d$SCOREstd <- scale(d$SCORE, center=TRUE, scale=TRUE)  
  colnames(d)[colnames(d)=="IID"] <- "id"
  colnames(d)[colnames(d)=="SCORE"] <- scorenam
  colnames(d)[colnames(d)=="SCOREstd"] <- scorestd
  colnames(d)[colnames(d)=="CNT"] <- cntnam
  colnames(d)[colnames(d)=="CNT2"] <- cnt2nam
  d1 = subset(d, select = -c(PHENO, FID) )
  assign(paste("S",i,sep=""),d1)
}

# drop d and d1 dataframes
rm(d,d1)

# create dataframe list
df.list <-list(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,outcomeData)

# merge the datafreames together on the common "id" column
data <- Reduce(function(x, y) merge(x, y, all=T, by="id"), df.list, accumulate=F)

# save the merged PRS and outcome dataset
write.table(data, file="pgc_scz_scores_with_outcomes.txt", row.names=F, col.names=T, quote=F, sep ="\t")  

