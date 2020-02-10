install.packages("qqman", repos='http://cran.us.r-project.org')
install.packages("GenABEL", repos='http://cran.us.r-project.org')

library(qqman, quiet=TRUE)
library(GenABEL, quiet=TRUE)

qqplotpval <- function(P, filename=NULL)
{
	require(GenABEL)
	l <- estlambda(P, method="median")
	nom <- paste("lambda = ", round(l$estimate, 3), sep="")
	if(!is.null(filename))
	{
		png(filename)
	}
	estlambda(P, method="median", plot=TRUE, main=nom)
	if(!is.null(filename))
	{
		dev.off()
	}
}

arguments <- commandArgs(T)
infile <- arguments[1]
outfile <- arguments[2]


cat("\n\nReading", infile, "\n")

a <- read.table(infile, he=T)
a$CHR <- as.numeric(a$CHR)
a <- subset(a, P != 0)



qqfilename <- paste(outfile, "_qqplot.png", sep="")
manhattanfilename <- paste(outfile, "_manhattan.png", sep="")

cat("Making", qqfilename, "\n")

qqplotpval(a$P, qqfilename)

cat("Making", manhattanfilename, "\n")

png(file=manhattanfilename)
manhattan(a, main=paste(basename(outfile), "Manhattan plot"))
dev.off()
