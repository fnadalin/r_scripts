# Created on 2020/02/11

# Perform DEA between two cells subsets defined by identity slots and identity
# values
# The same object can be provided in --object1 and --object2
# FIXME: if data integration was performed, the DEA is done on the transformed expression values!!!
# ---> FIXED: provide the assay explicitely

########################### DEFAULT PARAMETER VALUES ###########################

IDENT_SLOT1 <- "orig.ident"
IDENT_SLOT2 <- "orig.ident"
TEST <- "MAST"

LABEL1 <- "obj1"
LABEL2 <- "obj2"

MIN_PERC <- 0.1
LOGFC <- 0.25
LOGFC_FILT <- 0.5
ADJ_PVAL_FILT <- 0.05
MIN_PERC_FILT <- 0.1

TOP_DEG <- 50
GENES_EXCLUDE <- NULL

NO_FILTER <- FALSE
NO_HEATMAP <- FALSE
NO_VOLCANO <- FALSE
NO_SAVE <- FALSE


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--object1", type = "character",
        help="[REQUIRED] .Robj file containing the 1st Seurat object"),
    make_option("--object2", type = "character",
        help="[REQUIRED] .Robj file containing the 2nd Seurat object"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] directory where output will be stored"),
    make_option("--ident_slot1", type = "character", default = IDENT_SLOT1,
        help="ID slot in object1@meta.data"),
    make_option("--ident_slot2", type = "character", default = IDENT_SLOT2,
        help="ID slot in object2@meta.data"),
    make_option("--ident1", type = "character",
        help="[REQUIRED] comma-separated list of identities for condition 1"),
    make_option("--ident2", type = "character",
        help="[REQUIRED] comma-separated list of identities for condition 2"),
    make_option("--label1", type = "character", default = LABEL1,
        help="re-assign each selected cells from object1 this label"),	
    make_option("--label2", type = "character", default = LABEL2,
        help="re-assign each selected cells from object2 this label"),	
	make_option("--test", type = "character", default = TEST,
        help="comma-separated list of statistical tests for DEGs computation"),
	make_option("--min_perc", type = "double", default = MIN_PERC,
        help="minimum fraction of cells where the gene is detected in either one of the two sets [default=%default]"),
	make_option("--logFC", type = "double", default = LOGFC,
        help="minimum logFC between the two sets in order for a gene to be considered [default=%default]"),
	make_option("--logFC_filt", type = "double", default = LOGFC_FILT,
        help="minimum logFC between the two sets in order for a gene to be selected among the filtered DEGs [default=%default]"),
	make_option("--adj_pval_filt", type = "double", default = ADJ_PVAL_FILT,
        help="minimum adjusted p-value in order for a gene to be selected among the filtered DEGs [default=%default]"),
	make_option("--min_perc_filt", type = "double", default = MIN_PERC_FILT,
        help="minimum fraction of cells where the gene is detected to be selected among the filtered DEGs [default=%default]"),
	make_option("--genes_exclude", type = "character",
        help="file containing the list of genes to be excluded from the filtered DEGs"),
	make_option("--top_deg", type = "integer", default = TOP_DEG,
        help="number of top significant filtered DEGs to plot"),
	make_option("--no_filter", action = "store_true",
        help="do not filter the DEGs"),
	make_option("--no_heatmap", action = "store_true",
        help="do not create a heatmap plot"),
	make_option("--no_volcano", action = "store_true",
        help="do not create a volcano plot"),
	make_option("--no_obj_save", action = "store_true",
        help="do not save a new object")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$object1)) {
	write("Option --object is required\nTry --help for help", stderr()) 
	q()
} else {
	OBJECT1 <- opt$object1
}

if (is.null(opt$object2)) {
	write("Option --object is required\nTry --help for help", stderr()) 
	q()
} else {
	OBJECT2 <- opt$object2
}

if (is.null(opt$out_dir)) {
	write("Option --out_dir is required\nTry --help for help", stderr()) 
	q()
} else {
	OUT_DIR <- opt$out_dir
}

if (!is.null(opt$ident_slot1)) {
	IDENT_SLOT1 <- opt$ident_slot1
} 

if (!is.null(opt$ident_slot2)) {
	IDENT_SLOT2 <- opt$ident_slot2
} 

if (is.null(opt$ident1)) {
	write("Option --ident1 is required\nTry --help for help", stderr()) 
	q()	
} else {
	IDENT1 <- unlist(strsplit(opt$ident1, split=","))	
}

if (is.null(opt$ident2)) {
	write("Option --ident2 is required\nTry --help for help", stderr()) 
	q()	
} else {
	IDENT2 <- unlist(strsplit(opt$ident2, split=","))
}

if (!is.null(opt$label1)) {
	LABEL1 <- opt$label1
}

if (!is.null(opt$label2)) {
	LABEL2 <- opt$label2
}

if (!is.null(opt$test)) 
	TEST <- opt$test

if (!is.null(opt$min_perc)) 
	MIN_PERC <- opt$min_perc

if (!is.null(opt$logFC)) 
	LOGFC <- opt$logFC

if (!is.null(opt$logFC_filt)) 
	LOGFC_FILT <- opt$logFC_filt

if (!is.null(opt$adj_pval_filt)) 
	ADJ_PVAL_FILT <- opt$adj_pval_filt

if (!is.null(opt$min_perc_filt)) 
	MIN_PERC_FILT <- opt$min_perc_filt

if (!is.null(opt$genes_exclude)) 
	GENES_EXCLUDE <- c(as.matrix(read.table(opt$genes_exclude)))

if (!is.null(opt$top_deg)) 
	TOP_DEG <- opt$top_deg

if (!is.null(opt$no_filter))
	NO_FILTER <- TRUE

if (!is.null(opt$no_heatmap))
	NO_HEATMAP <- TRUE

if (!is.null(opt$no_volcano))
	NO_VOLCANO <- TRUE

if (!is.null(opt$no_obj_save))
	NO_SAVE <- TRUE


############################### EXPORT FUNCTIONS ###############################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)


######################### PRINT THE PARAMETERS TO FILE #########################

dir.create(OUT_DIR, showWarnings = FALSE)
sink(file = paste(OUT_DIR, "parameters.txt", sep = "/"), append = TRUE)
cat("============ Differential expression ============\n")
for (i in 1:length(opt))
	cat(paste(names(opt)[i], opt[[i]], "\n", sep="="))
sink()


################################## EXECUTION ###################################

object1 <- readRDS(OBJECT1)
Idents(object1) <- IDENT_SLOT1
object1_sel <- subset(x = object1, idents = IDENT1)
n1 <- length(colnames(x = object1_sel))

object2 <- readRDS(OBJECT2)
Idents(object2) <- IDENT_SLOT2
object2_sel <- subset(x = object2, idents = IDENT2)
n2 <- length(colnames(x = object2_sel))

object <- merge(x = object1_sel, y = object2_sel, merge.data = TRUE)
Idents(object) <- c(rep(LABEL1, n1), rep(LABEL2, n2))

if (!NO_SAVE)
	save(object, file=file.path(OUT_DIR, "object.Robj"))

DEA_OUT_DIR <- file.path(OUT_DIR, TEST)
dir.create(DEA_OUT_DIR, showWarnings = FALSE)

prefix <- paste(DEA_OUT_DIR, "/DEG_", LABEL1, "_", LABEL2, sep = "")
deg.table <- paste(prefix, "tsv", sep=".")
filt.table <- paste(prefix, "_filtered.tsv", sep="")
plot.name <- paste(prefix, "_filtered_heatmap.pdf", sep="")
volcano.plot.name <- paste(prefix, "_filtered_volcano.pdf", sep="")

GeneMarkersTableNEW(object = object, 
                    out.name = deg.table, 
                    ident.1 = LABEL1, 
                    ident.2 = LABEL2, 
                    test.use = TEST, 
                    min.pct = MIN_PERC, 
                    logFC = LOGFC)

if (!NO_FILTER) {
	FilterClusterGeneMarkers(input.table = deg.table, 
		                     output.table = filt.table, 
		                     logFC.filt = LOGFC_FILT,
		                     adjpval.filt = ADJ_PVAL_FILT,
		                     min.pct = MIN_PERC_FILT, 
		                     num = TOP_DEG, 
		                     genes.filter = GENES_EXCLUDE)
	if (!NO_HEATMAP) {
		object <- ScaleData(object, assay = "RNA")
		HeatmapPlotNEW(object = object, 
				       filt.table = filt.table, 
				       plot.name = plot.name,
				       width = 6)
	}

	if (!NO_VOLCANO) {
		title <- paste(LABEL1, LABEL2, sep=" vs ")
		VolcanoPlotFilter(all.table = deg.table,
				          filt.table = filt.table, 
				          plot.name = volcano.plot.name,
				          genes.filter = GENES_EXCLUDE,
				          title = title,
				          subtitle = paste(IDENT_SLOT1, IDENT_SLOT2, sep = "\n"))
	}
}

sessionInfo()


q()
	
