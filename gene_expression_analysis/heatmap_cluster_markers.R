
# ~/.R/R_3.6.3_no_x11/bin/R

DIR <- "analysis/perturb"
OBJECT <- file.path(DIR, "object_filt_NEW.rds")
TEST <- "MAST"

PCS <- 11
K <- 50
RES <- c(0.1, 0.2, 0.3)

library("Seurat")
library("ComplexHeatmap")
library("dplyr")
library("scales")
library("RColorBrewer")
library("circlize")

object <- readRDS(OBJECT)

for (res in RES) {

	CL <- paste0("clusters_pca", PCS, "_k", K, "_res", res)
	CL_ALL_DEG_DIR <- file.path(DIR, CL, TEST)

	Idents(object) <- object@meta.data[[CL]]
	cl_ident <- sort(unique(Idents(object)))
	cl_ident <- factor(cl_ident, levels = cl_ident)

	genes <- NULL
	for (cl in cl_ident) {
		deg_file <- file.path(CL_ALL_DEG_DIR, paste0("DEG_", TEST, "_cl", cl, "-all.tsv"))
		degs <- read.table(deg_file, sep = "\t", header = TRUE)
		g <- degs$geneID[degs$avg_log2FC > 0]
		g <- as.character(g[!(g %in% genes)])
		genes <- c(genes, g[1:10])
	}

	# retain only the scaled data for the genes that I'm actually using
	object <- ScaleData(object, features = genes)
	expr <- AverageExpression(object = object, features = genes, slot = "scale.data")
	M <- as.matrix(expr$RNA)
	for (i in 1:nrow(M)) {
		for (j in 1:ncol(M)) {
			if (M[i,j] > 2.5)
				M[i,j] <- 2.5
			if (M[i,j] < -2.5)
				M[i,j] <- -2.5
		}
	}
	colnames(M) <- cl_ident

	col_scaled <- colorRamp2(c(-1, 0, 2), cm.colors(3))
	col_cl <- hue_pal()(length(cl_ident))
	names(col_cl) <- cl_ident

	ha_column <- HeatmapAnnotation(df = data.frame(conditions = cl_ident), col = list(conditions = col_cl), show_annotation_name = FALSE, simple_anno_size = unit(0.2, "cm"), show_legend = FALSE)

	ht_scaled <- Heatmap(M, col = col_scaled, name = "Avg expr (z-score)", column_title = CL, top_annotation = ha_column,
		show_row_names = TRUE, row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        show_column_names = TRUE, column_names_gp = gpar(fontsize = 10), column_names_rot = 45,
		cluster_columns = FALSE, cluster_rows = FALSE, column_split = cl_ident, column_names_side = "top",
		show_heatmap_legend = TRUE, heatmap_legend_param = list(direction = "horizontal")) 

	pdf(file.path(CL_ALL_DEG_DIR, "heatmap_top10.pdf"), height = 1.3*(1+length(cl_ident)), width = 4.5)
	draw(ht_scaled, heatmap_legend_side = "bottom")
	dev.off()

}


sessionInfo()
q()

