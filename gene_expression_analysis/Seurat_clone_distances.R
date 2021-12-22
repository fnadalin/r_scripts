
# for a given sample, select the cells assigned to clones and compute pairwise distances across those cells
# distances are computed on the PC space (dr = dimensional reduction), where only the top PCs are considered (dim = # PCs),
# according to the automatic criterion for the estimation of stdev due to noise
# plot the distance density: intra-clone (bold line), inter-clone (thin lines), for each clone.

DISTANCE <- "euclidean"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <object.Rds> <dr> <dim> <sample> <out_dir>\n")
    q()
}

OBJECT <- args[1]
dr <- args[2]
dim <- args[3]
SAMPLE <- args[4]
OUTDIR <- args[5]

library("Seurat")
library("ComplexHeatmap")
library("scales") # for alpha()
library("RColorBrewer")
library("circlize")

# extract the clones for the current sample
object <- readRDS(OBJECT)
idx <- which(object@meta.data$sample.name == SAMPLE & object@meta.data$expr.GBC.list != "")
data <- as.matrix(object@reductions[[dr]]@cell.embeddings[idx,1:dim])
expr.GBC.list <- as.character(object@meta.data$expr.GBC.list[idx])
names(expr.GBC.list) <- colnames(object)[idx]

# compute distances between cells
distances <- dist(data, method = DISTANCE)
dist_matrix <- as.matrix(distances, diag = 0)
max_dist <- max(dist_matrix)

# compute the number of occurrences of each clone
clone_occurrence <- table(expr.GBC.list)
num_clones <- dim(clone_occurrence)
idx_clones <- order(clone_occurrence, decreasing = TRUE)

perc90 <- idx_clones[1:ceiling(num_clones/10)]
perc20  <- idx_clones[(num_clones-floor(num_clones/5)):num_clones]

plotfile <- file.path(OUTDIR, paste0("clone_distance_density_", SAMPLE, ".pdf"))
pdf(plotfile, width = 16, height = 2*ceiling(length(perc90)/4))
par(mfrow = c(ceiling(length(perc90)/4), 8))
for (i in perc90) {
    c1 <- names(clone_occurrence)[i]
    idx1 <- which(expr.GBC.list == c1)
    # collect the distances within the same clone
    dm <- dist_matrix[idx1,idx1]
    intra_clone_dist <- unlist(lapply(2:length(idx1), function(x) dm[x,1:(x-1)]))
    d_intra <- density(intra_clone_dist)
    sub <- paste("frequency =", clone_occurrence[i])
    # collect cross-clone distances for high frequency clones
    plot(0, 0, type = "n", main = c1, xlab = paste(DISTANCE, "distance"), ylab = "density", xlim = c(0,max_dist), ylim = c(0,0.3), sub = sub)
    for (j in perc90) {
        c2 <- names(clone_occurrence)[j]
        if (c1 != c2) {
            idx2 <- which(expr.GBC.list == c2)
            dm <- dist_matrix[idx1,idx2]
            inter_clone_dist <- c(dm)
            d_inter <- density(inter_clone_dist)
            lines(d_inter$x, d_inter$y, lwd = 1, col = alpha("red", alpha = 0.1))
        }
    }
    lines(d_intra$x, d_intra$y, lwd = 2)
    # collect cross-clone distances for low frequency clones
    plot(0, 0, type = "n", main = c1, xlab = paste(DISTANCE, "distance"), ylab = "density", xlim = c(0,max_dist), ylim = c(0,0.3), sub = sub)
    for (j in perc20) {
        c2 <- names(clone_occurrence)[j]
        if (c1 != c2) {
            idx2 <- which(expr.GBC.list == c2)
            dm <- dist_matrix[idx1,idx2]
            inter_clone_dist <- c(dm)
            d_inter <- density(inter_clone_dist)
            lines(d_inter$x, d_inter$y, lwd = 1, col = alpha("blue", alpha = 0.1))
        }
    }
    lines(d_intra$x, d_intra$y, lwd = 2)
}
dev.off()

q()

# sort the cells based on the number of occurrences of the clones they belong to
expr.GBC.list_occurrence <- clone_occurrence[match(expr.GBC.list, names(clone_occurrence))]
idx2 <- order(expr.GBC.list_occurrence, decreasing = TRUE)
dist_matrix_sorted <- as.matrix(dist_matrix[idx2,idx2])
dist_matrix_sorted_labels <- names(expr.GBC.list_occurrence)[idx2]
dist_matrix_sorted_labels <- factor(dist_matrix_sorted_labels, levels = unique(dist_matrix_sorted_labels))

# plot the top 1000 cells (in terms of clone occurrence) on a heatmap

col_scaled <- colorRamp2(c(0,10,40), c("blue","white","red"))
N <- 1000
# col_cl <- hue_pal()(length(unique(dist_matrix_sorted_labels[1:N])))
col_cl <- viridis_pal()(length(unique(dist_matrix_sorted_labels[1:N])))
names(col_cl) <- unique(dist_matrix_sorted_labels[1:N])

ha_col <- HeatmapAnnotation(df = data.frame(cells = dist_matrix_sorted_labels[1:N]), col = list(cells = col_cl), show_annotation_name = TRUE, 
                            simple_anno_size = unit(0.3, "cm"), show_legend = FALSE, which = "column")
ha_row <- HeatmapAnnotation(df = data.frame(cells = dist_matrix_sorted_labels[1:N]), col = list(cells = col_cl), show_annotation_name = TRUE, 
                            simple_anno_size = unit(0.3, "cm"), show_legend = FALSE, which = "row")
ht_scaled <- Heatmap(dist_matrix_sorted[1:N,1:N], col = col_scaled, name = paste(DISTANCE, "distance"), top_annotation = ha_col, left_annotation = ha_row,
	             show_row_names = FALSE, show_column_names = FALSE, 
		     cluster_columns = TRUE, cluster_rows = TRUE, show_row_dend = FALSE, show_column_dend = FALSE,
                     column_split = dist_matrix_sorted_labels[1:N], row_split = dist_matrix_sorted_labels[1:N],
                     cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                     row_title_rot = 0, column_title_rot = 90, row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
		     show_heatmap_legend = TRUE, heatmap_legend_param = list(direction = "horizontal")) 

pdf(paste0("heatmap_cell_distances_", dr, "_top", N, "_", SAMPLE, ".pdf"), height = 5, width = 5.5)
draw(ht_scaled, heatmap_legend_side = "bottom")
dev.off()

# compute the average distance matrix (mean = cells of the same clone)
n <- length(unique(dist_matrix_sorted_labels))
M_avg <- matrix(NA, nrow = n, ncol = n)
rownames(M_avg) <- colnames(M_avg) <- unique(dist_matrix_sorted_labels)
i <- 1
for (h in 1:n) {
    n1 <- expr.GBC.list_occurrence[idx2[i]]
    j <- 1
    for (k in 1:n) {
        n2 <- expr.GBC.list_occurrence[idx2[j]]
        M_avg[h,k] <- mean(dist_matrix_sorted[i:(i+n1-1),j:(j+n2-1)])
        j <- j + n2
    }
    i <- i + n1
}

# plot the average distances on a matrix

col_scaled <- colorRamp2(c(0,10,40), c("blue","white","red"))
col_cl <- viridis_pal()(length(unique(dist_matrix_sorted_labels)))
names(col_cl) <- unique(dist_matrix_sorted_labels)

ha_col <- HeatmapAnnotation(df = data.frame(clones = unique(dist_matrix_sorted_labels)), col = list(clones = col_cl), show_annotation_name = TRUE, 
                            simple_anno_size = unit(0.3, "cm"), show_legend = FALSE, which = "column")
ha_row <- HeatmapAnnotation(df = data.frame(clones = unique(dist_matrix_sorted_labels)), col = list(clones = col_cl), show_annotation_name = TRUE, 
                            simple_anno_size = unit(0.3, "cm"), show_legend = FALSE, which = "row")
ht_scaled <- Heatmap(M_avg, col = col_scaled, name = paste(DISTANCE, "distance"), top_annotation = ha_col, left_annotation = ha_row,
	             show_row_names = FALSE, show_column_names = FALSE, column_title = paste0(SAMPLE, " clones"),
		     cluster_columns = TRUE, cluster_rows = TRUE, show_row_dend = TRUE, show_column_dend = TRUE,
		     show_heatmap_legend = TRUE, heatmap_legend_param = list(direction = "horizontal")) 

pdf(paste0("heatmap_cell_distances_avg_", dr, "_", SAMPLE, ".pdf"), height = 5, width = 4.5)
draw(ht_scaled, heatmap_legend_side = "bottom")
dev.off()



sessionInfo()


q()


