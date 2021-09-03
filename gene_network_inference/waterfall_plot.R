
PATH="/data/users/fnadalin/curie/exp/DC_HIV1_scRNASeq/scripts/functions_Seurat_20180504"

OBJECT <- ""
id_slot <- "orig.ident"
cl_filter <- ""
renaming <- ""
palette <- c()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript scenic_waterfall_plot.R <scenic_folder> <out_folder> <object> [<id_slot> <cl_filter> <renaming>]\n\n")
	cat("<scenic_folder>     SCENIC folder (containing int/ and output/)\n")
	cat("<out_folder>        output folder where the plots are saved\n")
	cat("<object>            name of the file containing the Seurat object\n")
	cat("<id_slot>           identity slot for identifying cell groups (OPTIONAL)\n")
	cat("<cl_filter>         comma-separated list of cluster IDs to be filtered out (OPTIONAL)\n")
	cat("<renaming>          comma-separated IDs to be assigned to the SELECTED clusters, after filtering (OPTIONAL)\n")
	cat("                    Example: \"5,1,4,2,6,3\" means 0=>5, 1=>1, 2=>4, 3=>2, 4=>6, 5=>3\n\n")
	q()
}

# parse the input parameters
SCENIC_FOLDER <- args[1]
OUT_FOLDER <- args[2]
OBJECT <- args[3]
if (length(args) > 3) {
	id_slot <- args[4]
	if (length(args) > 4) {
		cl_filter <- args[5]	
		if (length(args) > 5) {
			renaming <- args[6]
		}
	}
}

# load the libraries
library("SCENIC")
library("RcisTarget")
library("scales")

source(paste(PATH, "/functions.R", sep=""))

dir.create(OUT_FOLDER, showWarnings = FALSE)
OUT_FOLDER <- normalizePath(OUT_FOLDER)

object <- eval(parse(text=load(OBJECT)))
new.cl.names <- unlist(strsplit(renaming, split=","))
new.cl.names <- new.cl.names[(1:length(new.cl.names)-1) %in% unique(sort(object@meta.data[[id_slot]]))]
renaming <- paste(new.cl.names, collapse=",")

if (renaming != "") object <- RenameClusters(object, renaming, id_slot, cl_filter)
ident <- object@meta.data[[id_slot]]
cluster.ids <- unique(sort(ident))

# select colors
palette <- hue_pal()(sum(!(cluster.ids %in% c("X"))))
names(palette) <- cluster.ids[!(cluster.ids %in% c("X"))]

# load the regulon/cell AUC matrix
setwd(SCENIC_FOLDER)
scenicOptions <- readRDS("int/scenicOptions.Rds") 
# aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC <- readRDS("int/3.4_regulonAUC.Rds")
auc <- getAUC(aucell_regulonAUC)

# remove cells not belonging to the selected clusters
auc <- auc[,!(ident %in% c("X"))]
ident <- ident[!(ident %in% c("X"))]

# extract the AUC threshold automatically selected by SCENIC
aucell_thresholds <- readRDS("int/3.5_AUCellThresholds.Rds")

cells <- colnames(auc)
regulons <- rownames(auc)
for (i in 1:length(regulons)) {
	df <- data.frame(cells = cells, AUC = auc[i,], cluster = ident)
	df <- df[order(df$AUC, decreasing=TRUE),]
	df <- df[order(df$cluster),]
	df$cells <- 1:length(cells)
	thr <- aucell_thresholds[[regulons[i]]]$aucThr$selected
	g <- ggplot(data = df, aes(x=cells, y=AUC)) + geom_bar(stat="identity", aes(color=cluster, fill=cluster)) + ggtitle(regulons[i]) + geom_hline(yintercept = thr, linetype = "dashed", color = "blue")
	out_pdf <- file.path(OUT_FOLDER, paste(gsub(" .*", "", regulons[i]), "pdf", sep="."))
	pdf(out_pdf, width=5, height=3)
	print(g)
	dev.off()
}

q()














