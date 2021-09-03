
# given a sequence A and a sequence B, compute:
# 1. the hamming distance between A and B
# 2. the hamming distance between A and rc(B)
# 3. the levenshtein distance between A and B
# 4. the levenshtein distance between A and rc(B)
# 5. the existence of a longest common substring of length > 4 between A and B
# 6. the existence of a longest common substring of length > 4 between A and rc(B)

SUBSTR_LENGTH <- 5
MIN_HAM_DIST <- MIN_LV_DIST <- 4

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <barcodes.txt> <out_dir>\n")
    cat("\n<barcodes.txt> one per line\n\n")
    q()
}

BARCODES <- args[1]
OUT_DIR <- args[2]

library("stringdist") # for seq_dist
library("Biostrings") # for reverseComplement
library("plyr") # for revalue

dir.create(OUT_DIR, showWarnings = FALSE)

bb <- drop(as.matrix(read.table(BARCODES)))
bb_rc <- unlist(lapply(bb, function(x) as.character(reverseComplement(DNAString(x)))))
bb_num <- unlist(lapply(bb, function(x) paste(revalue(unlist(strsplit(x, split = "")), c("A" = 1, "C" = 2, "G" = 3, "T" = 4)), collapse = "")))  
bb_rc_num <- unlist(lapply(bb_rc, function(x) paste(revalue(unlist(strsplit(x, split = "")), c("A" = 1, "C" = 2, "G" = 3, "T" = 4)), collapse = "")))  

M_ham <- M_ham_rc <- M_lv <- M_lv_rc <- matrix(1000, ncol = length(bb), nrow = length(bb))
M_sub_exists <- M_sub_rc_exists <- matrix(0, ncol = length(bb), nrow = length(bb))

rownames(M_ham) <- rownames(M_ham_rc) <- rownames(M_lv) <- rownames(M_lv_rc) <- rownames(M_sub_exists) <- rownames(M_sub_rc_exists) <- bb
colnames(M_ham) <- colnames(M_ham_rc) <- colnames(M_lv) <- colnames(M_lv_rc) <- colnames(M_sub_exists) <- colnames(M_sub_rc_exists) <- bb

for (i in 1:length(bb)) {
    v_i <- as.numeric(unlist(strsplit(bb_num[i], split = "")))
    for (j in 1:length(bb)) {
        v_j <- as.numeric(unlist(strsplit(bb_num[j], split = "")))
        v_rc_j <- as.numeric(unlist(strsplit(bb_rc_num[j], split = "")))
        M_ham[i,j] <- seq_dist(v_i, v_j, method = "hamming")
        M_ham_rc[i,j] <- seq_dist(v_i, v_rc_j, method = "hamming")
        M_lv[i,j] <- seq_dist(v_i, v_j, method = "lv")
        M_lv_rc[i,j] <- seq_dist(v_i, v_rc_j, method = "lv")
        subs <- 0
        subs_rc <- 0
        for (k in 1:(length(bb)-SUBSTR_LENGTH+1)) {
            sub_k <- paste(v_j[k:(k+SUBSTR_LENGTH-1)], collapse = "")
            sub_rc_k <- paste(v_rc_j[k:(k+SUBSTR_LENGTH-1)], collapse = "")
            subs <- subs + length(grep(sub_k, bb[i]))
            subs_rc <- subs_rc + length(grep(sub_rc_k, bb[i]))
        }
        M_sub_exists[i,j] <- as.numeric(subs > 0)
        M_sub_rc_exists[i,j] <- as.numeric(subs_rc > 0)
    }
}

write.table(M_ham, file = file.path(OUT_DIR, "M_ham.tsv"), quote = FALSE, sep = "\t")
write.table(M_ham_rc, file = file.path(OUT_DIR, "M_ham_rc.tsv"), quote = FALSE, sep = "\t")
write.table(M_lv, file = file.path(OUT_DIR, "M_lv.tsv"), quote = FALSE, sep = "\t")
write.table(M_lv_rc, file = file.path(OUT_DIR, "M_lv_rc.tsv"), quote = FALSE, sep = "\t")
write.table(M_sub_exists, file = file.path(OUT_DIR, paste0("M_sub", SUBSTR_LENGTH, "_exists.tsv")), quote = FALSE, sep = "\t")
write.table(M_sub_rc_exists, file = file.path(OUT_DIR, paste0("M_sub", SUBSTR_LENGTH, "_rc_exists.tsv")), quote = FALSE, sep = "\t")

M_good <- matrix(0, ncol = length(bb), nrow = length(bb))

for (i in 1:length(bb)) {
    for (j in 1:length(bb)) {
        M_good[i,j] <- (M_ham[i,j] >= MIN_HAM_DIST &
                        M_ham_rc[i,j] >= MIN_HAM_DIST &
                        M_lv[i,j] >= MIN_LV_DIST &
                        M_lv_rc[i,j] >= MIN_LV_DIST &
                        M_sub_exists[i,j] == 0 &
                        M_sub_rc_exists[i,j] == 0)
    }
}

M_good <- apply(M_good, 2, as.numeric)
rownames(M_good) <- colnames(M_good) <- bb

write.table(M_good, file = file.path(OUT_DIR, "M_GOOD.tsv"), quote = FALSE, sep = "\t")


q()



