THRESH <- 1 # propensity threshold difference for a clone to be considered "closer" to a lineage

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <graph_dir> <outdir> <k>\n")
    q()
}

GRAPH_DIR <- args[1]
OUTDIR <- args[2]
K <- as.numeric(args[3])

MATRIX <- file.path(GRAPH_DIR, "nn_graph.Rds")
LABELS <- file.path(GRAPH_DIR, "node_labels.txt")
LINEAGES <- file.path(OUTDIR, "lineages.tsv")

M <- readRDS(MATRIX)
clone_class <- read.table(LABELS, sep = "\t")[,1]

clone_count <- table(clone_class)
clone_count <- clone_count[order(clone_count, decreasing = TRUE)]
clone_list <- names(clone_count)

lineages_anchor <- read.table(LINEAGES, sep = "\t")
lineages_anchor <- lineages_anchor[!is.na(lineages_anchor[,2]),]
lineages <- rep(0,length(clone_list))
names(lineages) <- clone_list
lineages[match(lineages_anchor[,1],names(lineages))] <- lineages_anchor[,2]
lineages_list <- unique(lineages_anchor[,2])

# iteratively assign missing clones to lineages

N <- length(clone_class)
for (c in clone_list) {
    # compute the observed number of edges between clone c and each lineage
    C_obs <- rep(0,length(lineages_list))
    if (lineages[c] == 0) { # clone not yet assigned to lineage
        n <- 0
        y <- 0
        for (x in M@i) {
            n <- n+1
            if (clone_class[x+1] == c) {
               while (y+1 <= length(M@p) & M@p[y+1] < n) {
                    y <- y+1
                }
                if (y != x+1) { # exclude loops
                    c1 <- clone_class[y]
                    l <- lineages[c1] # compute propensity on original lineages
                    if (l != 0) {
                        idx <- which(lineages_list == l)
                        C_obs[idx] <- C_obs[idx]+1
                    }
                }
            }
        }
    }
    # compute the expected number of edges between clone c and each lineage
    C_exp <- c()
    n_i <- sum(clone_class == c)
    for (l in lineages_list) {
        n_j <- sum(clone_count[lineages == l])
        C_exp <- c(C_exp, n_i*K*n_j / N)
    } 
    # update the lineages by possibly including c
    propensity <- C_obs/C_exp
    names(propensity) <- lineages_list
    m1 <- max(propensity)
    idx1 <- which(propensity == m1)
    if (length(idx1) == 1) {
        m2 <- max(propensity[-idx])
        if (m2 < 1 & m1 > m2+THRESH) {
            # assign c to lineage l
            l <- lineages_list[idx1]
            lineages[c] <- l
        }
    } 
}
    
out <- file.path(OUTDIR, "lineages_augmented.tsv")
write.table(lineages, file = out, quote = FALSE, sep = "\t", col.names = FALSE)

q()


