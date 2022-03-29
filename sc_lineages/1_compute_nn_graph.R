args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <indir> <outdir> <sample_id> <mode> <k>\n")
    q()
}

INDIR <- args[1]
OUTDIR <- args[2]
SAMPLE <- args[3]
HVG <- args[4]
K <- as.numeric(args[5])

# select the PC space to be used to compute the nearest neighbors 
PCA <- paste0("pca_", HVG)
NPCS <- file.path(INDIR, HVG, "num_PCs.txt")
nPCs <- read.table(NPCS)[,1]

OBJECT <- file.path(INDIR, "object.Rds")

# load libraries

library("Seurat")

# subset the object to the sample and compute the nearest neighbors

object <- readRDS(OBJECT)
object <- subset(object, subset = sample.name == SAMPLE) 
object <- FindNeighbors(object, k.param = K, reduction = PCA, dims = 1:nPCs) 

M <- object@graphs$RNA_nn
class <- object@meta.data$expr.GBC.list

saveRDS(M, file = file.path(OUTDIR, "nn_graph.Rds"))
write.table(class, file = file.path(OUTDIR, "node_labels.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)


q()


