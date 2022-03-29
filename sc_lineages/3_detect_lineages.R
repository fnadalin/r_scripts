METHOD <- "ward.D2"
NCLUST <- 2

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <top_clones> <outdir> [<nclust>]\n")
    q()
}

TOP_CLONES <- args[1]
OUTDIR <- args[2]
if (length(args) > 3) {
    NCLUST <- as.numeric(args[3])
}

MATRIX <- file.path(OUTDIR, "pair_propensity_symm_filt.tsv")

M <- as.matrix(read.table(MATRIX, sep = "\t", header = TRUE))
clones <- rownames(M)

d <- dist(M)
hc <- hclust(d, method = METHOD)
split <- cutree(hc, k = NCLUST)

all_clones <- read.table(TOP_CLONES, sep = "\t")[,1]
all_split <- rep(NA, length(all_clones))
all_split[match(clones, all_clones)] <- split
df <- data.frame(all_clones, all_split)

CLUST <- file.path(OUTDIR, "lineages.tsv")
write.table(df, file = CLUST, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

q()

