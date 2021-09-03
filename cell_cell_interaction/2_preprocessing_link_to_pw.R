
# Take as input an underscore-separated list of protein partners and return a dash-separated list to be used for PW analysis

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <P1_P2> <pvalues.txt> <out.txt>\n")
    cat("\n<pop1_pop2>    input partner list with P1 expressed in pop1 and P2 expressed in pop2 (\"P1_P2\"), separated by underscore, one per line\n")
    cat("<pvalues.txt>  as returned by CellPhoneDB\n")
    cat("<out.txt>      out file containing the dash-separated lists\n\n")
    q()
}

INFILE <- args[1]
PVALUES <- args[2]
OUTFILE <- args[3]

pairs <- drop(as.matrix(read.table(INFILE, sep = "\t")))
pvalues <- read.table(PVALUES, sep = "\t", header = TRUE)

idx <- which(pvalues$interacting_pair %in% pairs)
M <- matrix("", nrow = length(idx), ncol = 2)
for (i in 1:length(idx)) {
    j <- idx[i]
    p1 <- as.character(pvalues$gene_a[j])
    p2 <- as.character(pvalues$gene_b[j])
    if (p1 == "") {
        p1 <- gsub("complex:", "", as.character(pvalues$partner_a[j]))
    }
    if (p2 == "") {
        p2 <- gsub("complex:", "", as.character(pvalues$partner_b[j]))
    }
    v <- c(p1, p2)
    M[i,] <- v
}

write.table(M, file = OUTFILE, sep = "-", quote = FALSE, col.names = FALSE, row.names = FALSE)

q()

