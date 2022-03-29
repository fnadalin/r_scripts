args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <clone_list> <graph_dir> <outdir> <k>\n")
    q()
}

CLONE_LIST <- args[1]
GRAPH_DIR <- args[2]
OUTDIR <- args[3]
K <- as.numeric(args[4])

MATRIX <- file.path(GRAPH_DIR, "nn_graph.Rds")
LABELS <- file.path(GRAPH_DIR, "node_labels.txt")
clone_list <- read.table(CLONE_LIST, sep = "\t")[,1]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

M <- readRDS(MATRIX)
class <- read.table(LABELS, sep = "\t")[,1]

idx <- which(clone_list %in% class)
clone_list <- clone_list[idx]

clone_idx <- match(class, clone_list, nomatch = 0)

# compute the number of direct edges between all possible clone pairs (non-symmetric matrix)

C_obs <- matrix(0, nrow = length(clone_list), ncol = length(clone_list))
colnames(C_obs) <- rownames(C_obs) <- clone_list

n <- 0
y <- 0
for (x in M@i) {
    n <- n+1
    i <- clone_idx[x+1]
    if (i > 0) {
        while (y+1 <= length(M@p) & M@p[y+1] < n) {
            y <- y+1
        }
        if (y != x+1) { # exclude loops
            j <- clone_idx[y]
            if (j > 0) {
                C_obs[i,j] <- C_obs[i,j]+1
            }
        }
    }
}

# compute the expected number of direct edges 

C_exp <- matrix(0, nrow = length(clone_list), ncol = length(clone_list))
colnames(C_exp) <- rownames(C_exp) <- clone_list

N <- length(class)
for (i in 1:nrow(C_exp)) {
    n_i <- sum(class == clone_list[i])
    for (j in 1:ncol(C_exp)) {
        n_j <- sum(class == clone_list[j])
#        C_exp[i,j] <- n_i*K*n_j / N
        n_j <- sum(class == clone_list[j]) - (i == j)
        C_exp[i,j] <- n_i*K*n_j / (N-1) # corrected version: exclude node i from the total count
    }
}

# compute the propensity

C_prop <- matrix(0, nrow = length(clone_list), ncol = length(clone_list))
colnames(C_prop) <- rownames(C_prop) <- clone_list

for (i in 1:nrow(C_exp)) {
    for (j in 1:ncol(C_exp)) {
        C_prop[i,j] <- C_obs[i,j] / C_exp[i,j]
    }
}

# symmetrize

C_prop_symm <- t(C_prop)*C_prop
C_prop_symm <- apply(C_prop_symm, 2, sqrt)

out <- file.path(OUTDIR, "observed_pairs.tsv")
write.table(C_obs, file = out, quote = FALSE, sep = "\t")
out <- file.path(OUTDIR, "expected_pairs.tsv")
write.table(C_exp, file = out, quote = FALSE, sep = "\t")
out <- file.path(OUTDIR, "pair_propensity.tsv")
write.table(C_prop, file = out, quote = FALSE, sep = "\t")
out <- file.path(OUTDIR, "pair_propensity_symm.tsv")
write.table(C_prop_symm, file = out, quote = FALSE, sep = "\t")

# filter out the clones that are not close to themselves (i.e. propensity < 1)

idx_filt <- which(diag(C_prop) >= 1)
C_prop_filt <- C_prop[idx_filt,idx_filt]
C_prop_symm_filt <- C_prop_symm[idx_filt,idx_filt]

out <- file.path(OUTDIR, "pair_propensity_filt.tsv")
write.table(C_prop_filt, file = out, quote = FALSE, sep = "\t")

out <- file.path(OUTDIR, "pair_propensity_symm_filt.tsv")
write.table(C_prop_symm_filt, file = out, quote = FALSE, sep = "\t")



q()


