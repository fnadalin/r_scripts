# Created on 03/04/2021

# 1. create custom gene sets
# 2. run GSEA 

RunGSEACustom <- function (table, out_prefix, collection, order = "d", title = "", adj_pval = 0.1, nPerm = 1000) {

	genes_list <- read.table(table, header=TRUE, sep="\t")

	####### Prepare the gene list ###########

	## feature 1: numeric vector
	geneList <- genes_list[,2]
	## feature 2: named vector
	names(geneList) <- genes_list[,1]
	## feature 3: decreasing order
	if (order == "a") 
		geneList <- 0-geneList
	geneList <- sort(geneList, decreasing = TRUE)
	geneList <- geneList[!is.na(names(geneList))]

	############### Run GSEA ################

	gsea <- GSEA(geneList = geneList, exponent = 1, nPerm = nPerm, minGSSize = 0, maxGSSize = INFTY, TERM2GENE = collection, pvalueCutoff = adj_pval, pAdjustMethod = "BH", verbose = FALSE)
	# gsea_res <- gsea@result[gsea@result$pvalue < pval,]

	write.table(gsea@result, file = paste(out_prefix, ".tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

	sign_gene_sets <- which(gsea@result$p.adjust < adj_pval)
	if (length(sign_gene_sets) == 0)
		return()


	######## Plot the enriched sets #########

	cat_to_show <- min(nrow(gsea@result), 15)
	height <- max(4,0.5+0.4*cat_to_show)
	title <- title
	if (title != "") 
		title <- paste(title, " - ", sep=" ")
	title <- paste(title, "GSEA (custom gene sets)", sep="")

	plot <- paste(out_prefix, ".pdf", sep="")
	pdf(plot, width=12, height=height)
	print(dotplot(gsea, showCategory=cat_to_show, title=title))
	dev.off()

	if (cat_to_show > 0) {
		for (i in 1:cat_to_show) {
			plot <- paste(out_prefix, "_" , gsea@result$ID[i], ".pdf", sep="")
			pdf(plot, width=7, height=7)
			print(gseaplot(gsea, geneSetID = gsea@result$ID[i], title=gsea@result$Description[i]))
			dev.off()
		}
	}

	return()
}

INFTY <- 100000000
ADJ_PVAL_FILT <- 1.1 # means no filter on p-value

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript RunCustomGSEA.R <geneSetFiles> <geneSetNames> <queryFile> <out_dir>\n")
	cat("\n<geneSetFiles>   list of files containing gene sets, one per line (each file should contain one gene per line)\n")
	cat("<geneSetNames>   list of gene sets IDs, one per line\n")
	cat("<queryFile>      table of genes to be scored with GSEA (col1 = gene, col2 = score, one per line)\n")
	cat("<out_dir>        where the output should be saved\n\n")
	q()
}

GENE_SET_FILES <- args[1]
GENE_SET_NAMES <- args[2]
queryFile <- args[3]
out_dir <- args[4]

library("GSEABase")
library("clusterProfiler")

dir.create(out_dir, showWarnings = FALSE)

geneSetFiles <- drop(as.matrix(read.table(GENE_SET_FILES)))
geneSetNames <- drop(as.matrix(read.table(GENE_SET_NAMES)))

if (length(geneSetFiles) != length(geneSetNames)) 
	q()

# create the gene sets

collection <- data.frame(ont = c(), gene = c())
for (i in 1:length(geneSetFiles)) {
	genes <- drop(as.matrix(read.table(geneSetFiles[i])))
	df <- data.frame(ont = rep(geneSetNames[i],length(genes)), gene = genes)
	collection <- rbind(collection, df)
}

# score the query list of genes

out_prefix <- file.path(out_dir, "GSEAoutput")
RunGSEACustom(table = queryFile, out_prefix = out_prefix, collection = collection, order = "d", title = "", adj_pval = 1, nPerm = 1000)


sessionInfo()
q()
