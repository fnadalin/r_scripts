
METADATA_SLOTS <- NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <object.Rds> <out_dir> <samples> <metadata> [<metadata_slots>]\n")
    cat("\n<object.Rds>     input Seurat object, with normalized counts\n")
    cat("<out_dir>        output directory where to store the normalized count and the metadata, to be provided as input to CellPhoneDb\n")
    cat("<samples>        comma-separated list of sample names (to be retrieved from object@meta.data$sample.name)\n")
    cat("<metadata>       metadata slot to extract the cell subsets from (should be present in object@meta.data)\n")
    cat("<metadata_slots> comma-separated list of slots in the metadata field to retain [OPTIONAL - default: all slots]\n\n")
    q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
SAMPLES <- unlist(strsplit(args[3], split = ","))
METADATA <- args[4]
if (length(args) > 4) {
    METADATA_SLOTS <- unlist(strsplit(args[5], split = ","))
}

library("Seurat")

object <- readRDS(OBJECT)
Idents(object) <- object@meta.data$sample.name
object <- subset(object, idents = SAMPLES)
if (!is.null(METADATA_SLOTS)) {
    Idents(object) <- object@meta.data[[METADATA]]
    object <- subset(object, idents = METADATA_SLOTS)
}

object <- NormalizeData(object, normalization.method = "RC") # normalize but do not log-transform
counts <- GetAssayData(object, slot = "data")
meta <- object@meta.data[[METADATA]]

df.counts <- cbind(Gene = rownames(counts), as.data.frame(counts))
df.meta <- data.frame(Cell = colnames(counts), cell_type = meta)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

write.table(df.counts, file = file.path(OUT_DIR, "counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df.meta, file = file.path(OUT_DIR, "meta.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


sessionInfo()
q()


