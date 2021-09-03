
########################### DEFAULT PARAMETER VALUES ###########################

ORG <- "Hs"
ORGNAME <- "human"
SPECIES <- "Homo sapiens"

TOP_DEG <- 100
LOG2FC <- 0.5
ADJ_PVAL <- 0.05
MIN_PERC <- 0.1
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
    make_option("--deg_table", type = "character",
        help="[REQUIRED] .tsv file containing the DEGs (with fields: geneID, avg_log2FC, p_val_adj, pct.1, pct.2)"),
    make_option("--out_prefix", type = "character",
        help="[REQUIRED] output prefix where the files will be written"),
    make_option("--min_perc", type = "double", default = MIN_PERC,
        help="minimum fraction of cells where the gene is detected in either one of the two sets [default=%default]"),
    make_option("--log2FC", type = "double", default = LOG2FC,
        help="minimum logFC between the two sets in order for a gene to be considered [default=%default]"),
    make_option("--adj_pval", type = "double", default = ADJ_PVAL,
        help="maximum adjusted p-value in order for a gene to be selected [default=%default]"),
    make_option("--top_deg", type = "integer", default = TOP_DEG,
        help="maximum number of DEGs to consider [default=%default]"),
    make_option("--ann_qval", type = "double", default = QVAL,
        help="maximum q-value to label a GO term or PW as significant [default=%default]"),
    make_option("--GO_level", type = "integer", default = LEVEL,
        help="maximum GO level to consider [default=%default]"),
    make_option("--GO_similarity", type = "double", default = SIMILARITY,
        help="minimum similarity level between two GO terms to consider them redundant [default=%default]"),
    make_option("--title", type = "character", default = TITLE,
        help="title to be added to the plot (to identify the TITLEarison) [default=%default]")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$deg_table)) {
    write("Option --deg_table is required\nTry --help for help", stderr()) 
    q()
} else {
    DEG_TABLE <- opt$deg_table
}

if (is.null(opt$out_prefix)) {
    write("Option --out_prefix is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_PREFIX <- opt$out_prefix
}

if (!is.null(opt$min_perc)) 
    MIN_PERC <- opt$min_perc

if (!is.null(opt$log2FC)) 
    LOG2FC <- opt$log2FC

if (!is.null(opt$adj_pval)) 
    ADJ_PVAL <- opt$adj_pval

if (!is.null(opt$top_deg)) 
    TOP_DEG <- opt$top_deg

if (!is.null(opt$ann_qval)) 
    QVAL <- opt$ann_qval

if (!is.null(opt$GO_level)) 
    LEVEL <- opt$GO_level

if (!is.null(opt$GO_similarity)) 
    SIMILARITY <- opt$GO_similarity

if (!is.null(opt$title)) 
    TITLE <- opt$title


######################### PRINT THE PARAMETERS TO FILE #########################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "../functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

OUT_DIR <- dirname(OUT_PREFIX)
dir.create(OUT_DIR, showWarnings = FALSE)
sink(file = file.path(OUT_DIR, "parameters.txt"), append = TRUE)
for (i in 1:length(opt))
    cat(paste(names(opt)[i], opt[[i]], "\n", sep="="))
sink()


################################## EXECUTION ###################################

SEL_DEGS <- paste0(OUT_PREFIX, "_input_genes.txt")

GOenrichPreprocessing(deg.table = DEG_TABLE, sel.degs = SEL_DEGS, logFC.filt = LOG2FC, num = TOP_DEG, org = "Hs", to.official = TRUE)

go_out <- paste0(OUT_PREFIX, "_GO")
reactome_out <- paste0(OUT_PREFIX, "_REACTOME")
kegg_out <- paste0(OUT_PREFIX, "_KEGG")
msd_out <- paste0(OUT_PREFIX, "_MSIGDB")

GOenrichNEW(table = SEL_DEGS, out_prefix = go_out, org = ORG, title = TITLE, types = c("BP"))
PWenrichNEW(table = SEL_DEGS, out_prefix = reactome_out, db = "reactome", org = ORG, orgname = ORGNAME, title = TITLE)
PWenrichNEW(table = SEL_DEGS, out_prefix = kegg_out, db = "kegg", org = ORG, orgname = ORGNAME, title = TITLE)
KEGGpostprocessingNEW(input.table = paste0(kegg_out, ".tsv"), output.table = paste0(kegg_out, "_conv.tsv"), org = ORG)
MSDenrich(table = SEL_DEGS, out_prefix = msd_out, org = ORG, species = SPECIES, title = TITLE, types = c("C6","C7"))
for (t in c("C6","C7")) {
	KEGGpostprocessingNEW(input.table = paste0(msd_out, "_", t, ".tsv"), output.table = paste0(msd_out, "_", t, "_conv.tsv"), org = ORG)	
}



sessionInfo()
q()





