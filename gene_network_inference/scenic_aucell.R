# Run gene network inference with SCENIC starting from a Seurat objects

ORG_ID <- "hgnc"

GENES_DETECTED_IN_CELLS <- 20

CORES <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript ScenicGeneNetwork.R <matrix_dir> <out_dir> <local_dir> [<cores>]\n\n")
	cat("<matrix_dir>               name of the folder containing \"colnames.txt\", \"matrix.mtx\", \"rownames.txt\"\n")
	cat("<out_dir>                  where SCENIC folders int/ and output/ are found\n")
	cat("<local_dir>                local directory inside <out_dir> where output int/ folder with AUCell results will be saved\n")
	cat(paste0("<cores>                    number of cores for AUCell [OPTIONAL - default: ", CORES, "]\n\n"))
	q()
}

# parse the input parameters
MATRIX_DIR <- args[1]
OUT_DIR <- args[2]
LOCAL_DIR <- args[3]
if (length(args) > 3) {
	CORES <- as.numeric(args[4])	
}

# load the libraries
library("Matrix")
library("SingleCellExperiment")
library("SCENIC")
library("RcisTarget")

# print the parameters to file
sink(file = paste(OUT_DIR, "parameters.txt", sep="/"))
cat(paste("OBJECT=", OBJECT, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("LOCAL_DIR=", LOCAL_DIR, "\n", sep=""))
cat(paste("GENES_DETECTED_IN_CELLS=", GENES_DETECTED_IN_CELLS, "\n", sep=""))
sink()

M <- readMM(file.path(MATRIX_DIR), "matrix.mtx")
colnames(M) <- drop(as.matrix(read.table(file.path(MATRIX_DIR), "colnames.txt")))
rownames(M) <- drop(as.matrix(read.table(file.path(MATRIX_DIR), "rownames.txt")))

setwd(OUT_DIR)
scenicOptions <- readRDS(file="int/scenicOptions.Rds")
scenicOptions@settings$org <- ORG_ID
scenicOptions@settings$dbDir <- DB_DIR
scenicOptions@settings$nCores <- CORES

# rename the files in aucell so that the output will not be overwritten
dir.create(file.path(LOCAL_DIR, "int"), showWarnings = FALSE, recursive = TRUE)
for (names in grep("aucell|genesKept", rownames(scenicOptions@fileNames$int), value = TRUE)) {
	filename <- scenicOptions@fileNames$int[names,]
	scenicOptions@fileNames$int[names,] <- file.path(LOCAL_DIR, filename)
}
saveRDS(scenicOptions, file=file.path(LOCAL_DIR, "int/scenicOptions.Rds")) 

# filter the expression matrix
exprMat <- as.matrix(object@raw.data)
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions, minSamples = GENES_DETECTED_IN_CELLS) # this step will also removes the genes that are not present in the databases

# Use AUCells to score the regulons on individual cells -> here retrieve all the cells
# 1. AUCell_buildRankings() to build the gene ranking for each cell
# 2. AUCell_calcAUC() to compute the AUC of each gene set (regulon) for each cell
#    the parameter aucMaxRank is set to k; the AUC for (c,r) is the fraction of genes in regulon r that are found in the top k genes for c,
#    where k is the 1st percentile in the distribution of the number of detected genes per cell
exprMat_log <- log2(exprMat+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipHeatmap = TRUE, skipTsne = TRUE) # plots are not needed


q()

