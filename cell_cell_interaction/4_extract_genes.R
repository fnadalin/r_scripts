# Take as input a list of protein partner names and the list of all proteins/complexes information (returned by CellPhoneDB)
# Return the list of genes associated to each complex

orglib <- "org.Hs.eg.db"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <P1_P2> <P2_P1> <CellPhone_outdir> <out_dir>\n")
    cat("\n<pop1_pop2>         file containing complex list with P1 expressed in pop1 and P2 expressed in pop2, as returned by CellPhone, one per line\n")
    cat("<pop2_pop1>         file containing complex list with P1 expressed in pop2 and P2 expressed in pop1, as returned by CellPhone, one per line\n")
    cat("<CellPhone_outdir>  output CellPhoneDb folder containing pvalues.txt and deconvoluted.txt\n\n")
    q()
}

SELECTED1 <- args[1]
SELECTED2 <- args[2]
CELLPHONE_DIR <- args[3]
OUT_DIR <- args[4]

sel1 <- drop(as.matrix(read.table(SELECTED1, sep = "\t")))
sel2 <- drop(as.matrix(read.table(SELECTED2, sep = "\t")))

pval <- read.table(file.path(CELLPHONE_DIR, "pvalues.txt"), sep = "\t", header = TRUE)
deco <- read.table(file.path(CELLPHONE_DIR, "deconvoluted.txt"), sep ="\t", header = TRUE)

sel <- c(sel1, sel2)
idx1 <- which(pval$interacting_pair %in% sel)
complex_name <- pval$interacting_pair[idx1]
complex_id <- pval$id_cp_interaction[idx1]
idx2 <- which(deco$id_cp_interaction %in% complex_id)
df <- data.frame(complex_name = complex_name[match(deco$id_cp_interaction[idx2],complex_id)], gene_name = deco$gene_name[idx2])

write.table(df, file = file.path(OUT_DIR, "complex_and_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(unique(as.character(deco$gene_name[idx2])), file = file.path(OUT_DIR, "genes_unique.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)


sessionInfo()
q()

