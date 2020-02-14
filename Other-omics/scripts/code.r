## eQTL and meQTL practical in R
## Genetic Epidemiology Short Course 
## March 2020
## Bristol Medical School, University of Bristol 

## Created by Matthew Suderman 
## matthew.suderman@bristol.ac.uk
## Objectives
## 1. Prepare the data for analysis.
## 2. Conduct an eQTL analysis.
## 3. Produce plots to visually inspect findings.
## 4. (Optional) Conduct an meQTL analysis.
## 5. (Optional) Identify trios of associated genes, CpG sites and SNPs.

## Preliminaries
## Login to BlueCrystal and initiate an interactive session.
## qsub -I -q teaching -l nodes=1:ppn=1,walltime=02:00:00
## Change to the directory set up for this practical:
## cd /panfs/panasas01/sscm/shortcourse/genetic.epidemiology/pract25_meqtl
## You will not be able to save files here, so we will create an output directory in your home directory:
## mkdir ~/pract25_output

## Data
## In this practical we will run eQTL and meQTL analyses on a small publicly available dataset:
## Wagner, et al. The relationship between DNA methylation, genetic and expression inter-individual variation in untransformed human fibroblasts. Genome Biol. 2014 Feb 20;15(2):R37. GEO link
## It includes: 
## * genotype (Illumina Human1M-Duov3 DNA Analysis BeadChip), 
## * DNA methylation (Illumina HumanMethylation450 BeadChip), and 
## * gene expression (Illumina HumanRef-8 v3.0 Expression Beadchip) profiles from 57 human skin fibroblast samples obtained from Coriell and McGill Cellbank.


## Analysis steps
## 1. Preprocessing.
## a. Population stratification.
## b. Gene expression outliers.
## c. Matching samples between datasets.
## d. Non-genetic variation removal.
## 2. Software and data.
## a. R package MatrixEQTL.
## b. Genetic data.
## c. Population stratification covariates.
## d. Expression data.
## e. Gene and SNP locations.
## 3. eQTL analysis.
## 4. Aftermath.
## a. Number of associations.
## b. Top associations.
## c. Association plotting.
## 5. (Optional) meQTL analysis.
## 6. (Optional) Associated trios.


## 1a. Preprocessing: Population stratification
## This is detected by visualizing genetic distances between all pairs of individuals using multi-dimensional scaling (MDS).
## We will use PLINK for these computations (http://pngu.mgh.harvard.edu/~purcell/plink/plink2.shtml). Please ensure that PLINK is installed before continuing.
## module add apps/plink-1.90
## Calculate genetic distances:
## plink --bfile data/snp --genome --out ~/pract25_output/plink
## The output plink.genome file in the current directory gives all pairwise distances (see the DST column).

## 1a. Preprocessing: Population stratification in 2-dimensions
## We represent these distances in two dimensions as follows:
## plink --bfile data/snp --read-genome ~/pract25_output/plink.genome --cluster --mds-plot 2 --out ~/pract25_output/plink
## This generates a file plink.mds that looks something like this:
##  FID       IID    SOL           C1           C2 
##    1   GM02704      0   0.00919443  -0.00202899 
##    2   GM02706      0   0.00508846 -0.000520485 
##    3   GM01650      0   0.00754305  0.000371449 
##    4   GM01653      0   0.00489205  -0.00105916 
##    5   GM02640      0   0.00624695 -0.000418025 
##    6   GM02641      0   0.00520069  -0.00199805 
## FID - Family ID; IID - Individual ID (or sample name); SOL - Assigned solution code; C1 - Position on first dimension; C2 - Position on second dimension.

## 1a. Preprocessing: Population stratification visually
## Start R.  First you'll need to load the module:   
## module add languages/R-3.2.0
## To start R, simply type 'R'.
## R
## We'll send outputs from R to the output directory we set up earlier.  For convenience, we'll save the location as a variable.
out.dir <- "~/pract25_output"
## We now plot the coordinates in plink.mds in R as follows: 
snp.mds <- read.table(file.path(out.dir, "plink.mds"), header=T, sep="", stringsAsFactors=F)
pdf(file.path(out.dir, "snp-mds.pdf"))
plot(snp.mds[,"C1"], snp.mds[,"C2"], main="SNP MDS plot", xlab="first", ylab="second", pch=19)
dev.off()
## Note:  You can retrieve your plot from the pract25_output folder.  It may appear a bit different.
## Question: Can you identify the two potential outlier samples in the plot?

## 1b. Preprocessing: Gene expression outliers
## We could use a similar approach to identify gene expression outliers. 
rna <- read.csv("data/rna-data.csv", row.names=1)
rna.dist <- 1-abs(cor(rna))
rna.mds <- cmdscale(rna.dist,eig=TRUE, k=2)
pdf(file.path(out.dir, "rna-mds.pdf"))
plot(rna.mds$points[,1], rna.mds$points[,2],
     main="RNA MDS plot", xlab="first", ylab="second", pch=19)
dev.off()

## 1c. Preprocessing: Matching samples between datasets
## All too often samples are mislabelled in the lab. Fortunately, it is possible to computationally detect mismatches between genetic and gene expression data using MixupMapper:
## http://genenetwork.nl/wordpress/mixupmapper/
## In this case, we will skip this step, assuming that the study authors haven't made such errors!
## It is possible to do the same for genetic and DNA methylation data, particularly for data generated using the Illumina HumanMethylation450 BeadChip. It actually measures a selection of SNPs in addition to the methylation levels.

## 1d. Preprocessing: Non-genetic variation removal
## Gene expression measurements can be influenced by factors other than genetic variation. To maximize power to detect associations between genes and genetic variants, we can remove all non-genetic variation from the gene expression data. To save time, we will not do that here. In case you're curious, here are the basic steps:
## 1. Calculate the principal components of the gene expression data.
## 2. Test associations between principal components and genetic variants.
## 3. Remove any components with strong genetic associations.
## 4. Remove variation of the remaining principal components from the gene expression data.

## 2a. Software and data: R package MatrixEQTL
## If MatrixEQTL has not already been installed, we install it as follows:
## install.packages("MatrixEQTL")
## We then load the library in R as follows:
library(MatrixEQTL)

## 2b. Software and data: Genetic data
## MatrixEQTL operates on SlicedData objects. These objects provide access to large datasets without using huge amounts of computer memory.
## The genetic data is contained in the file data/snp-data.csv. We load it as follows:
snp.data <- SlicedData$new()            ## create a new SlicedData object
snp.data$fileDelimiter <- ","           ## columns in the file are separated by commas
snp.data$fileOmitCharacters <- "NA"     ## missing values are represented by "NA"
snp.data$fileSkipRows <- 1              ## the first column provides row names
snp.data$fileSkipColumns <- 1           ## the first row provides column names
snp.data$fileSliceSize <- 2000          ## only load the data 2000 rows at a time to save memory
snp.data$LoadFile("data/snp-data.csv") 
## Loading the data takes about half a minute.

## 2c. Software and data: Population stratification covariates
## When we test associations between genetic variants, we will need to adjust for population stratification. This is encoded by the multi-dimensional scaling coordinates that we created earlier.
## As required by MatrixEQTL, we load them as SlicedData objects:
covariates.data <- SlicedData$new()
covariates.data$initialize(t(snp.mds[,c("C1","C2")]))

## 2d. Software and data: Expression data
## Next we convert the gene expression data to SlicedData objects and verify that samples match the genotype data samples.
rna.data <- SlicedData$new()            ## create a new SlicedData object
rna.data$fileDelimiter <- ","           ## columns in the file are separated by commas
rna.data$fileOmitCharacters <- "NA"     ## missing values are represented by "NA"
rna.data$fileSkipRows <- 1              ## the first column provides row names
rna.data$fileSkipColumns <- 1           ## the first row provides column names
rna.data$fileSliceSize <- 2000          ## only load the data 2000 rows at a time to save memory
rna.data$LoadFile("data/rna-data.csv") 
## Don't forget to make sure that the gene expression and genetic datasets have the samples in the same order!
## stopifnot(all(colnames(snp.data) == colnames(rna.data)))

## 2e. Software and data: Gene and SNP locations
## Question: If we tested associations between all SNP-gene pairs, how many tests would we have to perform?

## That's a lot of tests. To save time, we will only consider cis pairs (i.e. SNPs that are with 1 million bases of the gene). MatrixEQTL can do this for us if we provide SNP and gene locations in the genome:
snp.loc <- read.csv("data/snp-features.csv", row.names=1)
rna.loc <- read.csv("data/rna-features.csv", row.names=1)

## 3. eQTL analysis: p-value thresholds
## Question: MatrixEQTL requires a p-value threshold for significance. To adjust for the number of tests, we could divide 0.05 by the number of tests to be performed. Because we only test cis pairs, the exact number of difficult to estimate. Can you think of a simple lower bound and upper bound on the number of tests?

## We'll use the lower bound to estimate a p-value threshold so that we don't miss any associations with statistically significant p-values.
threshold <- 0.05/lower.bound

## 3. eQTL analysis: finally!
## Finally we are ready to run the eQTL analysis.
eqtl <- Matrix_eQTL_main(snps = snp.data,
                         gene = rna.data,
                         cvrt = covariates.data, ## population stratification
                         pvOutputThreshold = 0, ## consider only cis pairs
                         pvOutputThreshold.cis = threshold, 
                         output_file_name.cis = file.path(out.dir, "eqtl-cis.txt"),
                         snpspos = snp.loc, 
                         genepos = rna.loc[,c("geneid","chr","left","right")],
                         cisDist = 1e6, ## define cis as < 1Mb
                         useModel = modelLINEAR, ## test using linear models 
                         verbose = TRUE,
                         pvalue.hist = TRUE,
                         min.pv.by.genesnp = FALSE,
                         noFDRsaveMemory = FALSE)
## This takes about 30 seconds.
## 4a. Aftermath: Number of associations
## Here is the exact number of tests that were performed:
eqtl$cis$ntests
## [1] 4863929
## Question: Above we guessed 1 test per SNP. How many tests were performed on average for each SNP? What p-value threshold should we apply given this number?

## At this threshold (p < 1.03e-08), we identify hundreds of associations.
eqtls <- eqtl$cis$eqtls
eqtls <- eqtls[which(eqtls$pvalue < threshold),]
nrow(eqtls)
## [1] 468
## Question: How many eQTLs did we identify?

## 4b. Aftermath: Top associations
## We will plot the strongest association:
idx <- which.min(eqtls$pvalue)
eqtls[idx,]
##         snps         gene statistic       pvalue          FDR      beta
## 1 rs17014489 ILMN_2100085  21.33979 1.096262e-27 8.063919e-22 0.4323303

## We save the corresponding gene and SNP:
top.gene <- as.character(eqtls$gene[idx])
top.snp <- as.character(eqtls$snps[idx])
top.gene
## [1] "ILMN_2100085"
top.snp
## [1] "rs17014489"

## 4c. Aftermath: Association plotting (data extraction)
## To plot the association, we will need to extract the SNP and gene expression values for the pair. This can be a bit painful using SlicedData objects because they don't support direct access to individual rows or columns. Here we create our own function to do this. Don't worry about how it works, just copy and paste into R.
sliced.data.row <- function(object, row) {
    if (is.factor(row))
        row <- as.character(row)
    if (is.character(row)) {
        row <- which(rownames(object) == row)
        if (length(row) == 0)
            stop("Invalid row")
    }
    slice.size <- nrow(object[[1]])
    slice.idx <- floor(row/slice.size) + 1
    slice <- object[[slice.idx]]
    row <- row - slice.size * (slice.idx - 1)
    ret <- slice[row,]
    names(ret) <- colnames(object)
    ret             
}


## 4c. Aftermath: Association plotting
## We can use this function to obtain the genotypes for any given SNP. For example, below we extract data for the top SNP-gene pair.
genotypes <- sliced.data.row(snp.data, top.snp)
genotypes <- factor(genotypes, levels=0:2, labels=c("AA","AB","BB"))
expression.levels <- sliced.data.row(rna.data, top.gene)
## We plot their association as follows:
pdf(file.path(out.dir, "top-association.pdf"))
boxplot(expression.levels ~ genotypes,
        main=paste(top.snp,top.gene),
        ylim=range(expression.levels),
        outline=F)
stripchart(expression.levels ~ genotypes,
           method="jitter", add=TRUE, vertical=TRUE,
           col=c("blue","orange","black", pch=19))
dev.off()




## 5. (Optional) meQTL analysis
cpg.data <- SlicedData$new()
cpg.data$fileDelimiter <- ","           ## columns are separated by commas
cpg.data$fileOmitCharacters <- "NA"     ## missing values are represented by "NA"
cpg.data$fileSkipRows <- 1              ## the first column provides row names
cpg.data$fileSkipColumns <- 1           ## the first row provides column names
cpg.data$fileSliceSize <- 2000          ## only load the data 2000 rows at a time
cpg.data$LoadFile("data/cpg-data.csv") 

cpg.loc <- read.csv("data/cpg-features.csv", row.names=1)

threshold <- 0.05/nrow(snp.data)

meqtl <- Matrix_eQTL_main(snps = snp.data,
                          gene = cpg.data,
                          cvrt = covariates.data,
                          pvOutputThreshold = 0,  ## ignore trans associations
                          pvOutputThreshold.cis = threshold,
                          output_file_name.cis = file.path(out.dir, "meqtl-cis.txt"),
                          snpspos = snp.loc, 
                          genepos = cpg.loc[,c("geneid","chr","left","right")],
                          cisDist = 1e6,                  ## define cis as < 1Mb
                          useModel = modelLINEAR,
                          verbose = TRUE,
                          pvalue.hist = TRUE,
                          min.pv.by.genesnp = FALSE,
                          noFDRsaveMemory = FALSE)

threshold <- 0.05/meqtl$cis$ntests
meqtls <- meqtl$cis$eqtls
meqtls <- meqtls[which(meqtls$pvalue < threshold),]



## 6. (Optional) Associated trios
## Such trios are of interest because they could indicate cases where a SNP modifies methylation levels which in turn modify gene expression levels.
## Construct the trios:
trios <- merge(data.frame(snp=eqtls$snps, gene=eqtls$gene, eqtl.statistic=eqtls$statistic),
               data.frame(snp=meqtls$snps, cpg=meqtls$gene, meqtl.statistic=meqtls$statistic))
## There are 157 trios.
## Calculate p-values for each gene/CpG site association:
trios$p.value <- sapply(1:nrow(trios), function(i) {
    expr <- sliced.data.row(rna.data, trios$gene[i])
    meth <- sliced.data.row(cpg.data, trios$cpg[i])
    cor.test(meth, expr)$p.value
})

## Identify the strongest gene/CpG site association.
idx <- which.min(trios$p.value)
trios[idx,]
##          snp         gene eqtl.statistic        cpg meqtl.statistic
## 6 rs10220917 ILMN_1656045       8.606951 cg25118879       -12.20575
##        p.value
## 6 2.513491e-15
