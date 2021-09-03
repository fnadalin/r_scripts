
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript extract_sparse_submatrix_from_seurat_object.R <object> <cell_ids> <out_dir>\n\n")
	q()
}

# parse the input parameters
OBJECT <- args[1]
CELL_IDS <- args[2]
OUT_DIR <- args[3]

library("Seurat")
library("Matrix")

object <- readRDS(OBJECT)
cell_ids <- drop(as.matrix(read.table(CELL_IDS)))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

M <- GetAssayData(object, slot = "counts")
m <- M[,match(cell_ids, colnames(M))]

writeMM(m, file = file.path(OUT_DIR, "matrix.mtx"))
write.table(cell_ids, file = file.path(OUT_DIR, "colnames.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(rownames(M), file = file.path(OUT_DIR, "rownames.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

sessionInfo()
q()

