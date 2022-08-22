
###### SEURAT V3 #######

# 1. Select the top significant PCs
# 2. Run the clustering with different parameters (number of neighbors and resolution)

# TODO: run TSNE and UMAP on the top significant PCs

GetNPCs <- function(in_dir) {
        pc_file <- file.path(in_dir,"num_PCs.txt")
        nPCs <- drop(as.matrix(read.table(pc_file)))

	nPCs
}

DrPlot <- function(object, dr = "pca", dr_plot = "pca", pcs = 50, pval = 1e-05, non.rand.sd.frac = 0.5, k, res, out_dir, dr_type = "PCA") {

	pdf(file.path(out_dir, paste0(dr_type, "_plot_samples.pdf")), width = 3.7*2, height = 3*2)
	print(DimPlot(object, group.by = "sample.name", reduction = dr_plot))
	dev.off()
	
	pdf(file.path(out_dir, paste0(dr_type, "_plot_ccScore.pdf")), width = 3.7*2, height = 3*2)
	print(DimPlot(object, group.by = "Phase", reduction = dr_plot))
	dev.off()

	for (kk in k) {
		for (r in res) {
			cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
			pdf(file.path(out_dir, paste0(dr_type, "_plot_", cl_ident, ".pdf")), width = 3.7*2, height = 3*2)
			print(DimPlot(object, group.by = cl_ident, reduction = dr_plot))
			dev.off()
		}
	}

}



args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: <obj.Robj> <in_dir> <out_dir> <params>\n")
	cat("\n<obj.Robj>   input Seurat object created at step 0\n")
	cat("<in_dir>     input directory containing the list of genes (from previous step) and the results of PCs selection\n")
	cat("<out_dir>    output directory containing the plots\n")
	cat("<params>     file including the values for the parameters, separated by \"=\"\n\n")
	q()
}

OBJECT <- args[1]
IN_DIR <- args[2]
OUT_DIR <- args[3]
PARAMS <- args[4]


params <- as.matrix(read.table(PARAMS, sep="="))
cores <- as.numeric(params[params[,1]=="cores", 2])
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))

library("Seurat")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

dir.create(OUT_DIR, showWarnings = FALSE)

object <- readRDS(OBJECT)

for (feature_method in c("mean.var.plot", "vst")) {
	if (feature_method == "mean.var.plot") {
		for (ymin in disps) {
			case <- paste("mean.var.plot_disp", ymin, sep = "")
			dr <- paste("pca", case, sep = "_")
			tsne_dr <- paste0("tsne_", dr)
			umap_dr <- paste0("umap_", dr)
			in_dir <- file.path(IN_DIR, case)
			out_dir <- file.path(OUT_DIR, dr)
			dir.create(out_dir, showWarnings = FALSE)
			nPCs <- GetNPCs(in_dir)
			if (!(tsne_dr %in% names(object@reductions)))
				object <- RunTSNE(object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = paste0("tSNE_", dr, "_"))
			if (!(umap_dr %in% names(object@reductions)))
                                object <- RunUMAP(object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = paste0("UMAP_", dr, "_"))
			DrPlot(object, dr = dr, dr_plot = dr, k = k, res = res, out_dir = out_dir, dr_type = "PCA")	
			DrPlot(object, dr = dr, dr_plot = tsne_dr, k = k, res = res, out_dir = out_dir, dr_type = "TSNE")
                        DrPlot(object, dr = dr, dr_plot = umap_dr, k = k, res = res, out_dir = out_dir, dr_type = "UMAP")	
		}
	} else { 
		for (nfeat in nfeats) {
			case <- paste("vst_top", nfeat, sep = "")
			dr <- paste("pca", case, sep = "_")
			tsne_dr <- paste0("tsne_", dr)
			umap_dr <- paste0("umap_", dr)
			in_dir <- file.path(IN_DIR, case)
			out_dir <- file.path(OUT_DIR, dr)
			dir.create(out_dir, showWarnings = FALSE)
			nPCs <- GetNPCs(in_dir)
			if (!(tsne_dr %in% names(object@reductions)))
				object <- RunTSNE(object, reduction = dr, dims = 1:nPCs, reduction.name = tsne_dr, reduction.key = tsne_dr)
			if (!(umap_dr %in% names(object@reductions)))
                                object <- RunUMAP(object, reduction = dr, dims = 1:nPCs, reduction.name = umap_dr, reduction.key = umap_dr)
			DrPlot(object, dr = dr, dr_plot = dr, k = k, res = res, out_dir = out_dir, dr_type = "PCA")	
			DrPlot(object, dr = dr, dr_plot = tsne_dr, k = k, res = res, out_dir = out_dir, dr_type = "TSNE")
                        DrPlot(object, dr = dr, dr_plot = umap_dr, k = k, res = res, out_dir = out_dir, dr_type = "UMAP")
		}
	}
}

saveRDS(object, file=OBJECT)


sessionInfo()
q()

