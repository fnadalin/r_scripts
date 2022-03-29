# Created on 2020/05/28

# VERY VERY IMPORTANT!!!!
# FOR SOME REALLY WEIRD REASON, FINDMARKERS DOES NOT ACCEPT "MAST" AS SPECIFIED
# IN THE PREAMBLE OF THIS FILE, JUST BELOW, BUT ONLY IF IT IS EXPLICITLY PROVIDED
# WITH THE --test OPTION!!!!

########################### DEFAULT PARAMETER VALUES ###########################

NO_PAIRS <- FALSE

# TEST <- "MAST"

MIN_PERC <- 0.1
LOGFC <- 0.25
LOGFC_FILT <- 0.5
ADJ_PVAL_FILT <- 0.05

LATENT_VARS <- NULL

TOP_DEG <- 50
GENES_EXCLUDE <- NULL

CORES <- 8


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--object", type = "character",
        help="[REQUIRED] .Robj file containing the Seurat object"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] directory where output will be stored"),
    make_option("--cl_mode", type = "character",
        help="[REQUIRED] mode ID for clustering (stored in object@meta.data)"),
    make_option("--test", type = "character", 
        help="[REQUIRED] statistical test for DEGs computation"),
    make_option("--no_pairs", action = "store_true",
        help="do not run DEA between all possible cluster pairs"),
    make_option("--min_perc", type = "double", default = MIN_PERC,
        help="minimum fraction of cells where the gene is detected in either one of the two sets [default=%default]"),
    make_option("--logFC", type = "double", default = LOGFC,
        help="minimum logFC between the two sets in order for a gene to be considered [default=%default]"),
    make_option("--latent_vars", type = "character", default = LATENT_VARS,
        help="comma-separated latent variables to be included in the model"),    
    make_option("--logFC_filt", type = "double", default = LOGFC_FILT,
        help="minimum logFC between the two sets in order for a gene to be selected among the filtered DEGs [default=%default]"),
    make_option("--adj_pval_filt", type = "double", default = ADJ_PVAL_FILT,
        help="minimum adjusted p-value in order for a gene to be selected among the filtered DEGs [default=%default]"),
    make_option("--genes_exclude", type = "character",
        help="file containing the list of genes to be excluded from the filtered DEGs"),
    make_option("--top_deg", type = "integer", default = TOP_DEG,
        help="number of top significant filtered DEGs to plot")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$object)) {
    write("Option --object is required\nTry --help for help", stderr()) 
    q()
} else {
    OBJECT <- opt$object
}

if (is.null(opt$out_dir)) {
    write("Option --out_dir is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_DIR <- opt$out_dir
}

if (is.null(opt$cl_mode)) {
    write("Option --cl_mode is required\nTry --help for help", stderr()) 
    q()
} else {
    CL_MODE <- opt$cl_mode
}

if (is.null(opt$test)) {
    write("Option --test is required\nTry --help for help", stderr()) 
    q()
} else {
    TEST <- opt$test
}

NO_PAIRS <- (!is.null(opt$no_pairs))

if (!is.null(opt$min_perc)) 
    MIN_PERC <- opt$min_perc

if (!is.null(opt$logFC)) 
    LOGFC <- opt$logFC

if (!is.null(opt$latent_vars))
    LATENT_VARS <- unlist(strsplit(opt$latent_vars, split = ","))

if (!is.null(opt$logFC_filt)) 
    LOGFC_FILT <- opt$logFC_filt

if (!is.null(opt$adj_pval_filt)) 
    ADJ_PVAL_FILT <- opt$adj_pval_filt

if (!is.null(opt$genes_exclude)) 
    GENES_EXCLUDE <- c(as.matrix(read.table(opt$genes_exclude)))

if (!is.null(opt$top_deg)) 
    TOP_DEG <- opt$top_deg


######################### PRINT THE PARAMETERS TO FILE #########################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "../functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

dir.create(OUT_DIR, showWarnings = FALSE)
sink(file = file.path(OUT_DIR, "parameters.txt"), append = TRUE)
cat("============ Differential expression ============\n")
for (i in 1:length(opt))
    cat(paste(names(opt)[i], opt[[i]], "\n", sep="="))
sink()


################################## EXECUTION ###################################

DEA_OUT_DIR <- file.path(OUT_DIR, TEST)
dir.create(DEA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(OBJECT)

plan("multiprocess", workers = CORES)

genes <- rownames(GetAssayData(object, slot="counts"))[!(rownames(GetAssayData(object, slot="counts")) %in% GENES_EXCLUDE)]

if (!NO_PAIRS) {
    ClusterGeneMarkersByPairs(object = object, 
                              out.dir = DEA_OUT_DIR, 
                              id = CL_MODE, 
                              test.use = TEST, 
                              min.pct = MIN_PERC, 
                              logFC = LOGFC,
                              latent.vars = LATENT_VARS)
}
ClusterGeneMarkersVsAll(object = object, 
                        out.dir = DEA_OUT_DIR, 
                        id = CL_MODE, 
                        test.use = TEST, 
                        min.pct = MIN_PERC,
                        logFC = LOGFC,
                        latent.vars = LATENT_VARS)

sessionInfo()
q()

if (!NO_PAIRS) {
    FilterClusterGeneMarkersPairs(object = object, 
                                  out.dir = DEA_OUT_DIR, 
                                  id = CL_MODE, 
                                  test.use = TEST, 
                                  logFC.filt = LOGFC_FILT,
                                  adjpval.filt = ADJ_PVAL_FILT, 
                                  num = TOP_DEG,
                                  genes.use = genes,
                                  heatmap = TRUE,
                                  volcano = TRUE)
}
FilterClusterGeneMarkersAll(object = object, 
                            out.dir = DEA_OUT_DIR, 
                            id = CL_MODE, 
                            test.use = TEST,
                            logFC.filt = LOGFC_FILT, 
                            adjpval.filt = ADJ_PVAL_FILT, 
                            num = TOP_DEG,
                            genes.use = genes,
                            heatmap = TRUE,
                            volcano = TRUE)

sessionInfo()


q()


