
# Take as input a list of protein partners (das-separated as in preprocessing step) and the list of all proteins/complexes information (returned by CellPhoneDB)
# Find the PW that contain all the proteins for each partner pair (including the genes corresponding to protein complexes)

orglib <- "org.Hs.eg.db"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <P1_P2> <P2_P1> <deconvoluted.txt> <out_dir>\n")
    cat("\n<pop1_pop2>         file containing input partner list with P1 expressed in pop1 and P2 expressed in pop2 (\"P1-P2\"), one per line\n")
    cat("<pop2_pop1>         file containing input partner list with P1 expressed in pop2 and P2 expressed in pop1 (\"P2-P1\"), one per line\n")
    cat("<deconvoluted.txt>  list of all proteins analyzed, returned by CellPhoneDb, containing the elements of each pair\n")
    cat("<out_dir>           ouptut directory containing the table of partner pairs (rows) - PW (columns) associations, and PW names\n\n")
    q()
}

POP1_POP2 <- args[1]
POP2_POP1 <- args[2]
DECONVOLUTED <- args[3]
OUT_DIR <- args[4]

library("org.Hs.eg.db")
library("clusterProfiler")
library("reactome.db")

pop1_pop2 <- drop(as.matrix(read.table(POP1_POP2, sep = "\t")))
pop2_pop1 <- drop(as.matrix(read.table(POP2_POP1, sep = "\t")))
deconvoluted <- read.table(DECONVOLUTED, sep = "\t", header = TRUE)

# remove homodimer interactions
homodimers <- pop1_pop2[pop1_pop2 %in% pop2_pop1]
pop1_pop2 <- pop1_pop2[!(pop1_pop2 %in% homodimers)]
pop2_pop1 <- pop2_pop1[!(pop2_pop1 %in% homodimers)]

# initialize the lists of genes associated to each partner pair
gene_lists <- pw_lists <- list()
pop2_pop1_rev <- unlist(lapply(pop2_pop1, function(x) paste(rev(unlist(strsplit(x, split = "-"))), collapse = "-")))
names <- c(pop1_pop2, pop2_pop1_rev)

# scan the partner pairs
pairs <- c(pop1_pop2, pop2_pop1)
for (i in 1:length(pairs)) {
    v <- unlist(strsplit(pairs[i], split = "-"))
    # retreive gene IDs for every unit in the partners 
    geneID <- c()
    for (el in v) {
        if (sum(deconvoluted$gene_name == el) > 0) {
            geneID <- c(geneID, el)
        } else {
            idx <- which(deconvoluted$complex_name == el)
            els <- as.character(unique(deconvoluted$gene_name[idx]))
            geneID <- c(geneID, els)
        }
    }
    # convert to ENTREZ ID
    SymbToEnt <- bitr(geneID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orglib, drop = TRUE)
    entrez <- SymbToEnt$ENTREZID
    n_units <- length(entrez)
    # extract the PW that contain any of the genes encoding for the units of the current partners
    all_pw <- unlist(lapply(entrez, function(x) reactomeEXTID2PATHID[[x]]))
    unique_pw <- unique(all_pw)
    # retriece the PW that contain all the genes encoding for the units
    pw_occ <- unlist(lapply(unique_pw, function(x) sum(all_pw == x)))
    pw_all_genes <- unique_pw[pw_occ == n_units]
    # populate the lists
    gene_lists[[names[i]]] <- geneID
    pw_lists[[names[i]]] <- pw_all_genes
}

pw <- unique(unlist(pw_lists))
n_pw <- length(pw)

# put the info on a matrix
# M[i,j] = 1 if pair i have all genes in pw j; 0 otherwise
M <- matrix(0, nrow = length(pairs), ncol = n_pw)
for (i in 1:nrow(M)) {
    p <- pw_lists[[names[i]]]
    idx <- which(pw %in% p)
    M[i,idx] <- rep(1,length(idx))
}

rownames(M) <- pairs
colnames(M) <- pw
pw_names <- unlist(lapply(pw, function(x) gsub(".*: ", "", reactomePATHID2NAME[[x]])))
pw_size <- unlist(lapply(pw, function(x) length(reactomePATHID2EXTID[[x]])))

# print everything
OUT_TABLE <- file.path(OUT_DIR, "pairs_pw_matrix.tsv")
write.table(M, file = OUT_TABLE, quote = FALSE, sep = "\t")

OUT_PW_NAMES <- file.path(OUT_DIR, "pw_id_names.tsv")
df <- data.frame(pw, pw_names, pw_size)
write.table(df, file = OUT_PW_NAMES, quote = FALSE, row.names = FALSE, sep = "\t")


q()


