# Run gene network inference with SCENIC starting from a Seurat objects

DB_DIR <- Sys.getenv(x = "DB_DIR")
DB <- "hg19-500bp-upstream-7species.mc9nr.feather"
ORG_ID <- "hgnc"

PERPLEXITY <- 10

SEED <- 123
GENES_DETECTED_IN_CELLS <- 20 # a gene must be detected in at least these many cells
METHOD <- "top10perTarget"

CORES <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript ScenicGeneNetwork.R <matrix_dir> <out_dir> [<cores>]\n\n")
	cat("<matrix_dir>               name of the folder containing \"colnames.txt\", \"matrix.mtx\", \"rownames.txt\"\n")
	cat("<out_dir>                  where SCENIC folders int/ and output/ will be created\n")
	cat(paste0("<cores>                    number of cores [OPTIONAL - default: ", CORES, "]\n\n"))
	q()
}

# parse the input parameters
MATRIX_DIR <- args[1]
OUT_DIR <- args[2]
if (length(args) > 2) {
	CORES <- as.numeric(args[3])
}

# load the libraries
library("SingleCellExperiment")
library("Matrix")
library("SCENIC")
library("RcisTarget")
library("doRNG")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# print the parameters to file
sink(file = paste(OUT_DIR, "parameters.txt", sep="/"))
cat(paste("MATRIX_DIR=", MATRIX_DIR, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("CORES=", CORES, "\n", sep=""))
cat(paste("SEED=", SEED, "\n", sep=""))
cat(paste("GENES_DETECTED_IN_CELLS=", GENES_DETECTED_IN_CELLS, "\n", sep=""))
cat(paste("METHOD=", METHOD, "\n", sep=""))
sink()

M <- readMM(file.path(MATRIX_DIR, "matrix.mtx"))
colnames(M) <- drop(as.matrix(read.table(file.path(MATRIX_DIR, "colnames.txt"))))
rownames(M) <- drop(as.matrix(read.table(file.path(MATRIX_DIR, "rownames.txt"))))

# IMPORTANT!!!! Set the options after having moved to the current directory
# the directories int/ and output/ are created automatically
setwd(OUT_DIR)
myDatasetTitle <- paste("SCENIC analysis on", MATRIX_DIR) # choose a name for your analysis
scenicOptions <- initializeScenic(org = ORG_ID, dbDir = DB_DIR, db = DB, datasetTitle = myDatasetTitle, nCores = CORES)
scenicOptions@settings$seed <- SEED
scenicOptions@settings$defaultTsne$perpl <- PERPLEXITY
saveRDS(scenicOptions, file = "int/scenicOptions.Rds") 

# Now run SCENIC
# the intermediate files are automatically stored inside int/

# filter the expression matrix
exprMat <- as.matrix(M)
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions, minSamples = GENES_DETECTED_IN_CELLS) # this step will also removes the genes that are not present in the databases
exprMat_filtered <- exprMat[genesKept,]

# GENIE3
# the output is automatically stored inside int/
# WARNING!!!! THIS STEP TAKES A LOT OF TIME
set.seed(SEED)
runGenie3(exprMat_filtered, scenicOptions)

# the file with the results is stored into RDSname
# run the following commands to import the results
# RDSname <- getIntName(scenicOptions, "genie3ll")
# genie3 <- readRDS(RDSname)

# scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # select the DB to use; comment to use both
runCorrelation(exprMat_filtered, scenicOptions) # this function has been introduced with SCENIC v1.1.1-9

# create the network modules
# 1. only edges with IM > 0.005 (importance)
# 2. top 50 IM targets for each TF
# 3. top 5, 10, 50 TF for each target
# (see [Aibar Nat. Met. 2017])
# Here also create a binarized matrix that contains 1 for positive correlation (> 0.3), -1 for negative correlation (< -0.3), and 0 otherwise
runSCENIC_1_coexNetwork2modules(scenicOptions) 

# prune the network and create the regulons
# to choose themethod for selection:
# > tfModules_asDF <- loadInt(scenicOptions, "tfModules_asDF")
# > unique(tfModules_asDF$method)
# multiple methods can be used
# WARNING!!!! THIS STEP TAKES A LOT OF TIME
# Here only positive correlations are retained (see previous step)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=METHOD)

# Use AUCells to score the regulons on individual cells -> here retrieve all the cells
# 1. AUCell_buildRankings() to build the gene ranking for each cell
# 2. AUCell_calcAUC() to compute the AUC of each gene set (regulon) for each cell
#    the parameter aucMaxRank is set to k; the AUC for (c,r) is the fraction of genes in regulon r that are found in the top k genes for c,
#    where k is the 1st percentile in the distribution of the number of detected genes per cell
exprMat_log <- log2(exprMat+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipHeatmap = TRUE, skipTsne = TRUE) # do not plot the tSNE because of the error in X11 library loading


sessionInfo()
q()

