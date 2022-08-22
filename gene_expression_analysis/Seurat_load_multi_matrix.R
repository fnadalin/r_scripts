
# R_LIBS_USER="/home/francesca/.local/R-3.6.3/"
# export R_LIBS_USER

# Take a list of CellRanger output matrix folders as input and generate a Seurat object

MIN_GENES <- 0
MIN_CELLS <- 0
USE_CELL_LIST <- FALSE
VARS_TO_REGRESS <- NULL


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--matrix_dir_list", type = "character",
        help="[REQUIRED] tsv file containing the input matrix directories, one per line, where each column represents:\n		1 - the sample name\n		2 - the matrix directory (both CellRanger v2 and v3 are supported)\n		3 - the file with the list of cell IDs to retain, one per line [OPTIONAL]"),
    make_option("--out_RDS", type = "character",
        help="[REQUIRED] output seurat object containing the cell barcodes from the matrices"),
    make_option("--min_genes", type = "integer", default = MIN_GENES,
        help="[OPTIONAL] minimum number of detected genes to keep a cell [default=%default]"),
    make_option("--min_cells", type = "integer", default = MIN_CELLS,
        help="[OPTIONAL] minimum number of cells where a gene is detected to keep it [default=%default]"),
    make_option("--vars_to_regress", type = "character",
        help="[OPTIONAL] comma-separated list of variables to be regressed out upon scaling")
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

if (is.null(opt$out_RDS)) {
    write("Option --out_RDS is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_FILE <- opt$out_RDS
}

if (!is.null(opt$min_genes))
    MIN_GENES <- as.numeric(opt$min_genes)

if (!is.null(opt$min_cells))
    MIN_CELLS <- as.numeric(opt$min_cells)

if (!is.null(opt$vars_to_regress))
    VARS_TO_REGRESS <- unlist(strsplit(opt$vars_to_regress, split = ","))


################################## EXECUTION ###################################

library("Matrix")
library("Seurat")

matrix_dir_list <- as.matrix(read.table(MATRIX_DIR_LIST, sep = "\t"))
if (ncol(matrix_dir_list) < 2) {
    write("The file provided with --matrix_dir_list should contain at least 2 columns\nTry --help for help", stderr()) 
    q()
}

if (ncol(matrix_dir_list) > 2) 
    USE_CELL_LIST <- TRUE

meta_sample_names <- c()
meta_sample_id <- c()
for (i in 1:nrow(matrix_dir_list)) {

    SAMPLE_NAME <- matrix_dir_list[i,1]
    MATRIX_DIR <- matrix_dir_list[i,2]
    if (USE_CELL_LIST) 
        CELL_LIST <- matrix_dir_list[i,3]

    mtx <- file.path(MATRIX_DIR, "matrix.mtx")
    mtx_gz <- paste0(mtx,".gz")

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

    if (USE_CELL_LIST) {
        v <- drop(as.matrix(read.table(CELL_LIST, sep = "\t")))
        idx <- which(barcodes %in% v)
        barcodes <- barcodes[idx]
        M <- M[,idx]
    }
    cells <- paste(SAMPLE_NAME, barcodes, sep = "-")
    M <- M[features[,3] == "Gene Expression",]
    rownames(M) <- genes
    colnames(M) <- cells

    if (i == 1) {
        MM <- matrix(NA, nrow = length(genes), ncol = 0)
        rownames(MM) <- genes
        colnames(MM) <- c()
    }
    
    if (sum(genes %in% rownames(MM)) != sum(rownames(MM) %in% genes)) {
        write("The gene lists must be the same across the matrices!!!\nTry --help for help", stderr()) 
        q()
    }

    MM <- cbind(MM, M)
    meta_sample_names <- c(meta_sample_names, rep(SAMPLE_NAME, ncol(M)))
    meta_sample_id <- c(meta_sample_id, rep(i, ncol(M)))
}

meta <- data.frame(sample.id = meta_sample_id, sample.name = meta_sample_names)
rownames(meta) <- colnames(MM)

object <- CreateSeuratObject(counts = MM, meta.data = meta, min.cells = MIN_CELLS, min.features = MIN_GENES)
rm(MM)

object <- NormalizeData(object)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

object <- ScaleData(object, vars.to.regress = VARS_TO_REGRESS)

saveRDS(object, file = OUT_FILE)


sessionInfo()
q()


