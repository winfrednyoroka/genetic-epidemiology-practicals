library(SKAT)
Generate_SSD_SetID("APOC3_maf0.05.bed", "APOC3_maf0.05.bim", "APOC3_maf0.05.fam", "APOC3_maf0.05.SetID", 
"APOC3_maf0.05.SSD", "APOC3_maf0.05.Info")
FAM<-read.table("APOC3.cov",header=T,stringsAsFactors=F,na.strings="NA")

##############################################
# 	Open SSD, and Run SKAT
#
y<-FAM[,3]
# To use a SSD file, please open it first. 
# 
SSD.INFO<-Open_SSD("APOC3_maf0.05.SSD", "APOC3_maf0.05.Info")

# Number of samples 
SSD.INFO$nSample 

# Number of Sets
SSD.INFO$nSets

methods=c("davies","optimal.adj")

obj <- SKAT_Null_Model(y ~ 1, out_type = "C",n.Resampling=0, type.Resampling = "bootstrap")

SetID<-SSD.INFO$SetInfo$SetID
	Z<-Get_Genotypes_SSD(SSD.INFO, 1)
	
	pvalues <- SetID
	
	for (m in 1:length(methods))
	{
		re <- SKAT(Z, obj, kernel = "linear.weighted", method=methods[m], estimate_MAF=2)
		p.val1<-re$p.value
		pvalues <- rbind(pvalues,p.val1)	
	}


########

library(SKAT)
Generate_SSD_SetID("APOC3_maf0.01.bed", "APOC3_maf0.01.bim", "APOC3_maf0.01.fam", "APOC3_maf0.01.SetID", 
"APOC3_maf0.01.SSD", "APOC3_maf0.01.Info")
FAM<-read.table("APOC3.cov",header=T,stringsAsFactors=F,na.strings="NA")

##############################################
# 	Open SSD, and Run SKAT
#
y<-FAM[,3]
# To use a SSD file, please open it first. 
# 
SSD.INFO<-Open_SSD("APOC3_maf0.01.SSD", "APOC3_maf0.01.Info")

# Number of samples 
SSD.INFO$nSample 

# Number of Sets
SSD.INFO$nSets

methods=c("davies","optimal.adj")

obj <- SKAT_Null_Model(y ~ 1, out_type = "C",n.Resampling=0, type.Resampling = "bootstrap")

SetID<-SSD.INFO$SetInfo$SetID
	Z<-Get_Genotypes_SSD(SSD.INFO, 1)
	
	pvalues <- SetID
	
	for (m in 1:length(methods))
	{
		re <- SKAT(Z, obj, kernel = "linear.weighted", method=methods[m], estimate_MAF=2)
		p.val1<-re$p.value
		pvalues <- rbind(pvalues,p.val1)	
	}
	print(pvalues)
#####

	


library(SKAT)
Generate_SSD_SetID("APOC3_exonic.bed", "APOC3_exonic.bim", "APOC3_exonic.fam", "APOC3_exonic.SetID", 
"APOC3_exonic.SSD", "APOC3_exonic.Info")
FAM<-read.table("APOC3.cov",header=T,stringsAsFactors=F,na.strings="NA")

##############################################
# 	Open SSD, and Run SKAT
#
y<-as.numeric(FAM[,3])
# To use a SSD file, please open it first. 
# 
SSD.INFO<-Open_SSD("APOC3_exonic.SSD", "APOC3_exonic.Info")

# Number of samples 
SSD.INFO$nSample 

# Number of Sets
SSD.INFO$nSets

methods=c("davies","optimal.adj")

obj <- SKAT_Null_Model(y ~ 1, out_type = "C",n.Resampling=0, type.Resampling = "bootstrap")

SetID<-SSD.INFO$SetInfo$SetID
	Z<-Get_Genotypes_SSD(SSD.INFO, 1)
	
	pvalues <- SetID
		
	for (m in 1:length(methods))
	{
		re <- SKAT(Z, obj, kernel = "linear.weighted", method=methods[m], estimate_MAF=2)
		p.val1<-re$p.value
		pvalues <- rbind(pvalues,p.val1)	
	}
        
##        
	
###
library(SKAT)
Generate_SSD_SetID("APOC3_exonic.bed", "APOC3_exonic.bim", "APOC3_exonic.fam", "APOC3_exonic.SetID", 
"APOC3_exonic.SSD", "APOC3_exonic.Info")
FAM<-read.table("APOC3.cov",header=T,stringsAsFactors=F,na.strings="NA")

##############################################
# 	Open SSD, and Run SKAT
#
y<-as.numeric(FAM[,3])
# To use a SSD file, please open it first. 
# 
SSD.INFO<-Open_SSD("APOC3_exonic.SSD", "APOC3_exonic.Info")

# Number of samples 
SSD.INFO$nSample 

# Number of Sets
SSD.INFO$nSets

methods=c("davies","optimal.adj")
X<-as.numeric(FAM[,4])
obj <- SKAT_Null_Model(y ~ X, out_type = "C",n.Resampling=0, type.Resampling = "bootstrap")

SetID<-SSD.INFO$SetInfo$SetID
	Z<-Get_Genotypes_SSD(SSD.INFO, 1)
	
	pvalues <- SetID
	
	for (m in 1:length(methods))
	{
		re <- SKAT(Z, obj, kernel = "linear.weighted", method=methods[m], estimate_MAF=2)
		p.val1<-re$p.value
		pvalues <- rbind(pvalues,p.val1)	
	}
	
	
	
###



