
# R_LIBS_USER="/home/francesca/.local/R-3.6.3/"
# export R_LIBS_USER

# Take a list of CellRanger output matrix folders as input and generate a Seurat object


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--matrix_dir_list", type = "character",
        help="[REQUIRED] comma-separated list of input directories containing the matrix files, with genes as rows and cell barcodes as columns"),
    make_option("--sample_names", type = "character",
        help="[REQUIRED] comma-separated list of sample names, corresponding to the matrix list\n"),
    make_option("--out_RDS", type = "character",
        help="[REQUIRED] output seurat object containing the cell barcodes from the matrices")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$matrix_dir_list)) {
	write("Option --matrix_dir_list is required\nTry --help for help", stderr()) 
	q()
} else {
	MATRIX_DIR_LIST <- opt$matrix_dir_list
}

if (is.null(opt$sample_names)) {
	write("Option --sample_names is required\nTry --help for help", stderr()) 
	q()
} else {
	SAMPLE_NAMES <- opt$sample_names
}

if (is.null(opt$out_RDS)) {
	write("Option --out_RDS is required\nTry --help for help", stderr()) 
	q()
} else {
	OUT_FILE <- opt$out_RDS
}


################################## EXECUTION ###################################

library("Matrix")
library("Seurat")

matrix_dir_list <- unlist(strsplit(MATRIX_DIR_LIST, split = ","))
sample_names <- unlist(strsplit(SAMPLE_NAMES, split = ","))

meta_sample_names <- c()
meta_sample_id <- c()
for (i in 1:length(matrix_dir_list)) {

	MATRIX_DIR <- matrix_dir_list[i]
	mtx <- file.path(MATRIX_DIR, "matrix.mtx")
        mtx_gz <- file.path(MATRIX_DIR, "matrix.mtx.gz")

	if (!xor( file.exists(mtx), file.exists(mtx_gz) )) {
		write(paste("Either", mtx, "or", mtx_gz, "must exist"), stderr())
		q()
	} 

	compress <- !file.exists(mtx)

	if (!compress) {
		M <- readMM(mtx)
		features <- as.matrix(read.table(file.path(MATRIX_DIR, "genes.tsv"), sep="\t"))
		barcodes <- drop(as.matrix(read.table(file.path(MATRIX_DIR, "barcodes.tsv"))))
	} else { 
		M <- readMM(gzfile(paste0(mtx,".gz")))
		gz <- gzfile(file.path(MATRIX_DIR, "features.tsv.gz"))
		features <- as.matrix(read.table(gz, sep="\t"))
		gz <- gzfile(file.path(MATRIX_DIR, "barcodes.tsv.gz"))
		barcodes <- drop(as.matrix(read.table(gz)))
	}
	genes <- features[features[,3] == "Gene Expression",2]
	cells <- gsub("-.*", paste0("-",i), barcodes)
	M <- M[features[,3] == "Gene Expression",]
	rownames(M) <- genes
	colnames(M) <- cells

	if (i == 1) {
		MM <- matrix(NA, nrow = length(genes), ncol = 0)
		rownames(MM) <- genes
		colnames(MM) <- c()
	}
	
	if (genes != rownames(MM)) {
		write("The gene lists must be the same across the matrices!!!\nTry --help for help", stderr()) 
		q()
	}

	MM <- cbind(MM, M)
	meta_sample_names <- c(meta_sample_names, rep(sample_names[i], ncol(M)))
	meta_sample_id <- c(meta_sample_id, rep(i, ncol(M)))
}

object <- CreateSeuratObject(counts = MM)
meta <- data.frame(sample.id = meta_sample_id, sample.name = meta_sample_names)
rownames(meta) <- colnames(object)
object <- AddMetaData(object = object, metadata = meta)

saveRDS(object, file = OUT_FILE)


q()


