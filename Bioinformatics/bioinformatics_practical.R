##########################################################
### Bioinformatics practical
### July 2022
##########################################################

#install.packages(c("httr","jsonlite"))
#install.packages("enrichR")
#install.packages("epigraphdb")

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#BiocManager::install("biomaRt")


library(httr)
library(jsonlite)
library(GenomicRanges)
library(biomaRt)
library(enrichR)
library(epigraphdb)



# DATASET
# Pairo-Castineira, E., Clohisey, S., Klaric, L. et al. Genetic mechanisms of critical illness in COVID-19.
# Nature 591, 92–98 (2021). https://doi.org/10.1038/s41586-020-03065-y
covid = matrix( c( "rs73064425", 3, 45901089, "LZTFL1",
		"rs9380142", 6, 29798794, "HLA-G",
		"rs143334143", 6, 31121426, "CCHCR1",
		"rs3131294", 6, 32180146, "NOTCH4",
		"rs10735079", 12, 113380008, "OAS1–OAS3",
		"rs2109069", 19, 4719443, "DPP9",
		"rs74956615", 19, 10427721, "TYK2",
		"rs2236757", 21, 34624917, "IFNAR2",
		"rs71325088", 3, 45862952, "LZTFL1",
		"rs6489867", 12, 113363550, "OAS1–OAS3",
		"rs11085727", 19, 10466123, "TYK2",
		"rs13050728", 21, 34615210, "IFNAR2" ),
		ncol=4, byrow=TRUE )
colnames(covid) = c("SNP","chr","pos","gene_symbol")
covid = as.data.frame(covid)



##########################################################
### PART 1 ###############################################
##########################################################


###
### 1.1 Look up your variants in dbSNP
###
#https://www.ncbi.nlm.nih.gov/snp/rs73064425
#https://www.ncbi.nlm.nih.gov/snp/rs9380142


###
### 1.2.a Load a reference genome and map the variants to gene.
###

mapping_nearest = c("nearest_gene","nearest_gene_symbol","nearest_distance")

#retrieve all genes with their GRCh37 coordinates from biomart
mart_grch37 = useEnsembl(biomart="ensembl",GRCh=37)
mart_grch37 = useDataset("hsapiens_gene_ensembl", mart_grch37)

attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","hgnc_symbol")
filters = c("chromosome_name")
values = list(chromosome_name=c(as.character(1:22),'X','Y') )

grch37 = getBM(attributes=attributes, filters=filters, values=values, mart=mart_grch37)

#this is our genome of reference
grch37 = GRanges(paste0("chr",grch37$chromosome_name),
		IRanges(as.numeric(grch37$start_position), as.numeric(grch37$end_position)), NULL,
		ensembl=grch37$ensembl_gene_id, symbol=grch37$hgnc_symbol )

grch37


# map our SNPs to the nearest gene in our reference
for (i in 1:nrow(covid)) {

	gr = GRanges(paste0("chr",covid$chr[i]), IRanges(as.numeric(covid$pos[i]), as.numeric(covid$pos[i])),  )
	nearest_gene = grch37[nearest(gr, grch37)]
	nearest_gene_distance = distanceToNearest(gr, grch37)@elementMetadata@listData$distance

	mapping_nearest = rbind(mapping_nearest, c(nearest_gene$ensembl, nearest_gene$symbol, nearest_gene_distance) )
}


colnames(mapping_nearest) = mapping_nearest[1,]
mapping_nearest = as.data.frame(mapping_nearest[-1,])
rownames(mapping_nearest) = covid$SNP

mapping_nearest


###
### 1.2.b Have all mappings been successful? Manually correct the ones that failed.
###

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=ENSG00000257452
mapping_nearest["rs10735079","nearest_gene_symbol"] = "OAS1"
mapping_nearest["rs6489867","nearest_gene_symbol"] = "OAS1"

# AP000295.1 or IFNAR2 ?
#http://asia.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000159110;r=21:33229901-33265675
mapping_nearest["rs2236757","nearest_gene"] = "ENSG00000159110"
mapping_nearest["rs2236757","nearest_gene_symbol"] = "IFNAR2"

mapping_nearest


###
### 1.2.c What does the nearest_distance column tell us?
###


###
### 1.2.d How does our list of genes compare with the one shown in the paper?
###
cbind(SNP=covid$SNP, covid=covid$gene_symbol, nearest=mapping_nearest$nearest_gene_symbol)
#https://www.genecards.org/cgi-bin/carddisp.pl?gene=RAVER1
#https://www.genecards.org/cgi-bin/carddisp.pl?gene=TYK2


###
### 1.3.a Request the eQTL associated to our variants.
###
eqtl = NULL
	
for (i in 1:nrow(covid)) {

	mapping_tmp = c("eqtl_gene","eqtl_tissue","eqtl_pval")
	
	cat("\rSNP: ", covid$SNP[i], " (", round(100*i/nrow(covid)), "%) ", sep="")

	request = GET(url = "http://www.ebi.ac.uk/eqtl/api/associations", 
		      query = list(
		              variant_id = covid$SNP[i],
		              size = 1000,
		              p_upper = 1E-8)
		     )
        
	if (request$status_code==200) {
	# 200 Ok

		# repeat till we get all eQTLs (API limits 1000 results), then we will select the one with lowest pval
		repeat{ 
			response = httr::content(request, as = "text", encoding = "UTF-8")
			df = jsonlite::fromJSON(response, flatten = TRUE)$`_embedded`$associations

			for (idx in names(df)) {
				mapping_tmp = rbind(mapping_tmp, as.character(df[[idx]][c("gene_id", "qtl_group", "pvalue")]) )
			}
			
			nextq = jsonlite::fromJSON(response, flatten = TRUE)$`_links`
			if( any(names(nextq)=="next") ) {
				request = GET(nextq$`next`$href)
				if (request$status_code!=200) break
			} else {
				break
			}

		}

		# is there any eQTL?
		if (!is.null(nrow(mapping_tmp)) ) {

			colnames(mapping_tmp) = mapping_tmp[1,]
		
			# sort by pval
			mapping_tmp = mapping_tmp[order(as.numeric(mapping_tmp[,3])),]

			eqtl[[covid$SNP[i]]] = mapping_tmp
		}
	}
}


###
### 1.3.b How many genes are associated to each variant? Are these associations present in all tissues?
###
eqtl[['rs9380142']][,c("eqtl_gene","eqtl_tissue")]




###
### 1.3.c Map to each SNP the gene with stronger association in lung, or if not present in blood.
###
mapping_eqtl = c("eqtl_gene", "eqtl_tissue", "eqtl_pval", "eqtl_gene_symbol")


for (i in 1:nrow(covid)) {

	if ( !is.null(eqtl[[covid$SNP[i]]]) ) {
	
		# eQTL in lung (associations are sorted by pval)
		sel = which(eqtl[[covid$SNP[i]]][,"eqtl_tissue"] == "Lung")[1]
		if (!is.na(sel)) {
			mapping_eqtl = rbind(mapping_eqtl, c(eqtl[[covid$SNP[i]]][sel,], NULL) )
			
		} else {
		
			# eQTL in blood (associations are sorted by pval)
			sel = which(eqtl[[covid$SNP[i]]][,"eqtl_tissue"] == "blood")[1]
			if (!is.na(sel)) {
				mapping_eqtl = rbind(mapping_eqtl, c(eqtl[[covid$SNP[i]]][sel,], NULL) )
			} else {
				mapping_eqtl = rbind(mapping_eqtl, NA )
			}
		}
				
	} else {
		mapping_eqtl = rbind(mapping_eqtl, NA )
	}
	
}

colnames(mapping_eqtl) = mapping_eqtl[1,]
mapping_eqtl = as.data.frame(mapping_eqtl[-1,])
rownames(mapping_eqtl) = covid$SNP


# retrieve gene symbols using biomart (eQTL Catalog returns ensembl)
attributes = c("ensembl_gene_id","hgnc_symbol")
filters = c("ensembl_gene_id")
values = list(ensembl_gene_id=mapping_eqtl$eqtl_gene)
mart_query = getBM(attributes=attributes, filters=filters, values=values, mart=mart_grch37)

for (i in 1:nrow(mapping_eqtl)) {

	sel = which( mart_query$ensembl_gene_id==mapping_eqtl$eqtl_gene[i] )[1]

	mapping_eqtl$eqtl_gene_symbol[i] = mart_query$hgnc_symbol[sel]

}

mapping_eqtl


###
### 1.3.d Check the mappings and, if needed, manually correct the ones that failed.
###

# https://www.genecards.org/cgi-bin/carddisp.pl?gene=ENSG00000272501
mapping_eqtl["rs143334143","eqtl_gene_symbol"] = "AL662844.4"

mapping_eqtl

###
### 1.3.c How does the eQTL genes compare with the previous mappings?
###

cbind(SNP=covid$SNP, covid=covid$gene_symbol, nearest=mapping_nearest$nearest_gene_symbol, eQTL=mapping_eqtl$eqtl_gene_symbol)




##########################################################
### PART 2 ###############################################
##########################################################


###
### 2.1 Look up at gene function
###
#https://www.genecards.org/cgi-bin/carddisp.pl?gene=LZTFL1#function
#https://www.genecards.org/cgi-bin/carddisp.pl?gene=HLA-G#function 


###
### 2.2.a What gene set libraries available in Enricher?
### For further analysis please select Reactome pathways and KEGG for human and mouse. 
###

listEnrichrDbs()

dbs = c("Reactome_2016", "KEGG_2019_Human", "KEGG_2019_Mouse")



###
### 2.2.b Using the nearest gene mapping, what pathways are enriched for our gene set?
###

# remove NA and duplicate genes
list_genes = mapping_nearest$nearest_gene_symbol
list_genes = unique(list_genes[!is.na(list_genes)])

enriched = enrichr(list_genes, dbs)

# for each database print enriched terms
for (j in 1:length(dbs)) {

	cat("\n",dbs[j],"\n", sep="")	
			
	sel = which(enriched[[j]]$Adjusted.P.value<0.05)
	for (i in sel) {
		p_name = enriched[[j]]$Term[i]
		p_genes = unlist(strsplit(enriched[[j]]$Genes[i], ";"))

		if (length(p_genes)>2) {
			cat("[",enriched[[j]]$Adjusted.P.value[i],"] ",p_name, ": ", sep="")	
			cat(p_genes, "\n", sep=" ")	
		}
	}

}




###
### 2.2.c Repeat the same analysis for the the eQTL mapping.
###

# remove NA and duplicate genes
list_genes = mapping_eqtl$eqtl_gene_symbol
list_genes = unique(list_genes[!is.na(list_genes)])

enriched = enrichr(list_genes, dbs)

# for each database print enriched terms
for (j in 1:length(dbs)) {

	cat("\n",dbs[j],"\n", sep="")	
			
	sel = which(enriched[[j]]$Adjusted.P.value<0.05)
	for (i in sel) {
		p_name = enriched[[j]]$Term[i]
		p_genes = unlist(strsplit(enriched[[j]]$Genes[i], ";"))

		if (length(p_genes)>2) {
			cat("[",enriched[[j]]$Adjusted.P.value[i],"] ",p_name, ": ", sep="")	
			cat(p_genes, "\n", sep=" ")	
		}
	}

}




##########################################################
### PART 3 ###############################################
##########################################################


###
### 3.1.a Find alternative gene targets that interact with our candidate gene HLA-DRB5.
###

endpoint = "/gene/druggability/ppi"
params = list(gene_name = "HLA-DRB5")
ppi_df = query_epigraphdb(route = endpoint, params = params, mode = "table")

ppi_df


###
### 3.1.b Filter those candidate genes with approved or clinical trial-phase drug candidates.
###

# Tier 1 genes encode protein targets of approved or clinical trial-phase drug candidates
sel = which(ppi_df$g2.druggability_tier == "Tier 1")
gene_list = c("HLA-DRB5",ppi_df$g2.name[sel])

gene_list


###
### 3.2.a Find evidence in the literature if our candidate genes are found to be associated with lung disease.
###

endpoint = "/gene/literature"
literature_df=NULL
# search for literature evidence for each candidate gene
for (i in 1:length(gene_list)) {
	params = list(
		gene_name = gene_list[i],
		object_name = "lung" # phenotype/disease
	    )
	if (is.null(literature_df))  literature_df = query_epigraphdb(route = endpoint, params = params, mode = "table")
	else  literature_df = rbind(literature_df, query_epigraphdb(route = endpoint, params = params, mode = "table"))
	
}

# expand the list of pubmed ids, refomartting
ref_ids = literature_df$pubmed_id
for (i in 1:nrow(literature_df)) {
	ids = toString(unlist(ref_ids[[i]]), collpase=",")
	literature_df$pubmed_id[i] = ids
}
literature_df$pubmed_id=unlist(literature_df$pubmed_id)

r = data.frame( gene = literature_df$gene.name,
		predicate = literature_df$st.predicate,
		disease = literature_df$lt.name,
		pubmed_id = literature_df$pubmed_id)

r


###
### 3.2.a Repeat the same analysis looking for associations with influenza.
###

endpoint = "/gene/literature"
literature_df=NULL
# search for literature evidence for each candidate gene
for (i in 1:length(gene_list)) {
	params = list(
		gene_name = gene_list[i],
		object_name = "influenza" # phenotype/disease
	    )
	if (is.null(literature_df))  literature_df = query_epigraphdb(route = endpoint, params = params, mode = "table")
	else  literature_df = rbind(literature_df, query_epigraphdb(route = endpoint, params = params, mode = "table"))
	
}


# expand the list of pubmed ids, refomartting
ref_ids = literature_df$pubmed_id
for (i in 1:nrow(literature_df)) {
	ids = toString(unlist(ref_ids[[i]]), collpase=",")
	literature_df$pubmed_id[i] = ids
}
literature_df$pubmed_id=unlist(literature_df$pubmed_id)

r = data.frame( gene = literature_df$gene.name,
		predicate = literature_df$st.predicate,
		disease = literature_df$lt.name,
		pubmed_id = literature_df$pubmed_id)

r




##########################################################
### END ##################################################
##########################################################

sessionInfo()


