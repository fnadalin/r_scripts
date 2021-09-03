
########################### DEFAULT PARAMETER VALUES ###########################

ORG <- "Hs"
ORGNAME <- "human"
SPECIES <- "Homo sapiens"

QVAL <- 0.05
LEVEL <- 4
SIMILARITY <- 0.7
TITLE <- "annotation"


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--genes", type = "character",
        help="[REQUIRED] .txt file containing the gene names, one per line"),
    make_option("--out_prefix", type = "character",
        help="[REQUIRED] output prefix where the files will be written"),
    make_option("--ann_qval", type = "double", default = QVAL,
        help="maximum q-value to label a GO term or PW as significant [default=%default]"),
    make_option("--GO_level", type = "integer", default = LEVEL,
        help="maximum GO level to consider [default=%default]"),
    make_option("--GO_similarity", type = "double", default = SIMILARITY,
        help="minimum similarity level between two GO terms to consider them redundant [default=%default]"),
    make_option("--title", type = "character", default = TITLE,
        help="title to be added to the plot (to identify the comparison) [default=%default]")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$genes)) {
    write("Option --genes is required\nTry --help for help", stderr()) 
    q()
} else {
    GENES <- opt$genes
}

if (is.null(opt$out_prefix)) {
    write("Option --out_prefix is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_PREFIX <- opt$out_prefix
}

if (!is.null(opt$ann_qval)) 
    QVAL <- opt$ann_qval

if (!is.null(opt$GO_level)) 
    LEVEL <- opt$GO_level

if (!is.null(opt$GO_similarity)) 
    SIMILARITY <- opt$GO_similarity


######################### PRINT THE PARAMETERS TO FILE #########################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

OUT_DIR <- dirname(OUT_PREFIX)
dir.create(OUT_DIR, showWarnings = FALSE)
sink(file = file.path(OUT_DIR, "parameters.txt"), append = TRUE)
for (i in 1:length(opt))
    cat(paste(names(opt)[i], opt[[i]], "\n", sep="="))
sink()


################################## EXECUTION ###################################

SEL_GENES <- paste0(OUT_PREFIX, "_input_genes.txt")

genes <- drop(as.matrix(read.table(GENES)))
orglib <- paste0("org.", ORG, ".eg.db")

# convert to official gene symbol
geneID <- alias2SymbolTable(genes, species = ORG) # to convert gene aliases to official gene names

# convert to Entrex ID
# NB: redundant terms are collapsed, this includes NA as well
SymbToEnt <- bitr(geneID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orglib, drop = FALSE)

# create the df with the 3 gene nomenclatures
df <- data.frame(Entrez = rep(NA,length(geneID)), geneID = geneID, geneIDorig = genes, group = rep(NA, length(geneID)))
df$Entrez[match(SymbToEnt$SYMBOL, df$geneID)] <- SymbToEnt$ENTREZID
write.table(df, file = SEL_GENES, sep="\t", quote = FALSE, row.names = FALSE)

go_out <- paste0(OUT_PREFIX, "_GO")
reactome_out <- paste0(OUT_PREFIX, "_REACTOME")
kegg_out <- paste0(OUT_PREFIX, "_KEGG")
msd_out <- paste0(OUT_PREFIX, "_MSIGDB")

GOenrichNEW(table = SEL_GENES, out_prefix = go_out, org = ORG, title = TITLE, types = c("BP"))
PWenrichNEW(table = SEL_GENES, out_prefix = reactome_out, db = "reactome", org = ORG, orgname = ORGNAME, title = TITLE)
PWenrichNEW(table = SEL_GENES, out_prefix = kegg_out, db = "kegg", org = ORG, orgname = ORGNAME, title = TITLE)
KEGGpostprocessingNEW(input.table = paste0(kegg_out, ".tsv"), output.table = paste0(kegg_out, "_conv.tsv"), org = ORG)
MSDenrich(table = SEL_GENES, out_prefix = msd_out, org = ORG, species = SPECIES, title = TITLE, types = c("C6","C7"))
for (t in c("C6","C7")) {
	KEGGpostprocessingNEW(input.table = paste0(msd_out, "_", t, ".tsv"), output.table = paste0(msd_out, "_", t, "_conv.tsv"), org = ORG)	
}



sessionInfo()
q()





