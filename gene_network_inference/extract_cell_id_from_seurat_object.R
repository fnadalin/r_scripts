
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript extract_cell_id_from_seurat_object.R <object> <out_file>\n\n")
	q()
}

# parse the input parameters
OBJECT <- args[1]
OUT_FILE <- args[2]

library("Seurat")

object <- readRDS(OBJECT)
write.table(colnames(object), file = OUT_FILE, quote = FALSE, col.names = FALSE, row.names = FALSE)

sessionInfo()
q()

