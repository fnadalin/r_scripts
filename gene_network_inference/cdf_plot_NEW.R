
# generate plots of the fraction of cells with AUC > x 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript cdf_plot.R <df_aucell_rds> <out_folder> <df_cell_info> <id_slot>\n\n")
	cat("<df_aucell_rds>     RDS object containing the regulon / cell matrix with AUC values\n")
	cat("<out_folder>        output pdf file where the plots are saved\n")
	cat("<df_cell_info>      tsv file containing cell information: cluster.id, sample.id... [MUST CONTAIN the \"cell.id\" field]\n")
	cat("<id_slot>           slot for identifying cell groups (must match <df_cell_info> header)\n\n")
	q()
}

# parse the input parameters
AUC_INFO <- args[1]
OUT_PDF <- args[2]
CELL_INFO <- args[3]
id_slot <- args[4]

library("scales")

# load the cell ID slot
cell_info <- read.table(CELL_INFO, sep = "\t", header = TRUE)

# load the regulon/cell AUC matrix
auc <- readRDS(AUC_INFO)
cell_info <- cell_info[cell_info$cell.id %in% colnames(auc),]
auc <- auc[order(rownames(auc)),match(cell_info$cell.id,colnames(auc))] # NEW: match the cell ID

# NEW: select only the cells that are present in the auc
ident <- cell_info[[id_slot]]
cluster.ids <- unique(sort(ident))
cluster.ids 

# select colors 
palette <- hue_pal()(length(cluster.ids))
names(palette) <- cluster.ids

regulons <- rownames(auc)
n <- length(regulons)

pdf(OUT_PDF, width = 8.27, height = 11.7)
par(mfrow=c(5,4))
plot_legend <- TRUE
for (i in 1:n) {
	df <- data.frame(AUC = drop(as.matrix(auc[i,])), cluster = ident)
	auc_steps <- min(df$AUC)+(0:100-0.5)*(max(df$AUC)-min(df$AUC))/100 # NEW: consider negative values as well
	cls <- sort(unique(df$cluster))
	plot(x = c(min(df$AUC),max(df$AUC)), y=c(0,1), pch = NA, xlab = "AUC", ylab = "fraction of cells", main = regulons[i]) # NEW: consider negative values as well
	for (j in 1:length(cls)) {
		AUC <- df$AUC[df$cluster == cls[j]]
		auc_cdf <- ecdf(AUC)
		lines(auc_steps, 1-auc_cdf(auc_steps), lwd = 1.5, col = palette[j]) # NEW: consider negative values as well
	}
	if (plot_legend)
		legend("topright", legend = cls, col = palette, lwd = 1.5, bty = "n")
	plot_legend <- FALSE
}
dev.off()

q()



