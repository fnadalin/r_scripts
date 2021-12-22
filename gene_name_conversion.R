
ORGLIB <- "org.Hs.eg.db"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <gene_list.txt> <from> <to> <out.txt>\n")
    cat("\n<gene_list.txt>    input file with gene names, one per line\n")
    cat("<from>,<to>        can be: ACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY, ONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT\n\n")
    q()
}

GENE_LIST <- args[1]
FROM <- args[2]
TO <- args[3]
OUT <- args[4]

library("org.Hs.eg.db")
library("clusterProfiler")

gene_list <- drop(as.matrix(read.table(GENE_LIST)))
converted <- bitr(gene_list, fromType = FROM, toType = TO, OrgDb = ORGLIB, drop=FALSE)
write.table(converted, file = OUT, quote = FALSE, sep = "\t", row.names = FALSE)

sessionInfo()
q()



