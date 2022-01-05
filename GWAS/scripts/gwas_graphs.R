# 04.01.22 Updated by R. Granell to remove GenABEL
# Problem reading input file: need to modify heading row 
# Problem using png command


 
install.packages("qqman", repos='http://cran.us.r-project.org')
#install.packages("GenABEL", repos='http://cran.us.r-project.org')

library(qqman, quiet=TRUE)
#library(GenABEL, quiet=TRUE)

#qqplotpval <- function(P, filename=NULL)
#{
#	require(GenABEL)
#	l <- estlambda(P, method="median")
#	nom <- paste("lambda = ", round(l$estimate, 3), sep="")
#	if(!is.null(filename))
#	{
#		png(filename)
#	}
#	estlambda(P, method="median", plot=TRUE, main=nom)
#	if(!is.null(filename))
#	{
#		dev.off()
#	}
#}

arguments <- commandArgs(T)
infile <- arguments[1]
outfile <- arguments[2]

qqfilename <- paste(outfile, "_qqplot.png", sep="")
manhattanfilename <- paste(outfile, "_manhattan.png", sep="")

cat("\n\nReading", infile, "\n")

# plink2 glm has produced an
# input file with  a heading row that starts with #. we need to remove this 
# also name of chromosome column is now CHROM
a <- read.table(infile, header=T, stringsAsFactors=F)
a<-na.omit(a)
a$CHR <- as.numeric(a$CHR)
a <- subset(a, P != 0)


qqfilename <- paste(outfile, "_qqplot.png", sep="")
manhattanfilename <- paste(outfile, "_manhattan.png", sep="")


#Calculate lambda
pvalue <- a$P
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda

cat("Making", qqfilename, "\n")
#QQ-plot
png(file=qqfilename)
#Add lambda in the tittle
#qq(a$P, main=paste(basename(outfile), "QQ-plot lambda=",lambda))
qq(a$P, main=paste(basename(outfile), "QQ-plot lambda="))
#qq(a$P)
dev.off()
exit
cat("Making", manhattanfilename, "\n")
#Manhattan plot
png(file=manhattanfilename)
manhattan(a,snp= "RSID", chr="CHR", bp="BP", p="p", main=paste(basename(outfile), "Manhattan plot"))
dev.off()




