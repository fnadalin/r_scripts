
# 1. Select features with different parameters (method and dispersion threshold / number of top features)
# 2. Compute the PCA 
# 3. Calculate the p-value of each PC 


Dr <- function(object, genes.list, dr = "pca", pcs = 50) {

	key <- paste(gsub('\\.|_',"",dr), "_", sep="")
	write("RunPCA", file="")
	object <- RunPCA(object, features = genes.list, npcs = pcs, verbose = FALSE, reduction.name = dr, reduction.key = key)

	# JackStraw only supports "pca"
#	names(object@reductions)[names(object@reductions) == dr] <- "pca"
#	write("JackStraw", file="")
#	object <- JackStraw(object, dims = pcs, num.replicate = 100)
#	write("ScoreJackStraw", file="")
#	object <- ScoreJackStraw(object, dims = 1:pcs)
#	names(object@reductions)[names(object@reductions) == "pca"] <- dr

	return(object)
}



args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: <obj.rds> <out_dir> <params>\n")
	cat("\n<obj.rds>    input Seurat object containing normalized and scaled data\n")
	cat("<out_dir>    output directory containing the results of feature selection and PCA\n")
	cat("<params>     file including the values for the parameters, separated by \"=\"\n\n")
	q()
}

library("Seurat")
library("future")

OBJECT <- args[1]
OUT_DIR <- args[2]
PARAMS <- args[3]


params <- as.matrix(read.table(PARAMS, sep="="))
cores <- as.numeric(params[params[,1]=="cores", 2])
xmin <- as.numeric(params[params[,1]=="xmin", 2])
xmax <- as.numeric(params[params[,1]=="xmax", 2])
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))


plan("multiprocess", workers = cores)

dir.create(OUT_DIR, showWarnings = FALSE)

object <- readRDS(OBJECT)

for (feature_method in c("mean.var.plot", "vst")) {
	if (feature_method == "mean.var.plot") {
		for (ymin in disps) {
			case <- paste0("mean.var.plot_disp", ymin)
			write(case, file="")
			dr <- paste("pca", case, sep="_")
			object <- FindVariableFeatures(object, selection.method = "mean.var.plot", mean.cutoff = c(xmin, xmax), dispersion.cutoff = c(ymin, Inf))
			out_dir <- file.path(OUT_DIR, case)
			dir.create(path = out_dir, showWarnings = FALSE)
			write.table(VariableFeatures(object), file = file.path(out_dir, "features.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
			object <- Dr(object, genes.list = VariableFeatures(object), dr = dr)
		}
	} else { 
		for (nfeat in nfeats) {
			case <- paste0("vst_top", nfeat)
			write(case, file="")
			dr <- paste("pca", case, sep="_")
			object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeat)
			out_dir <- file.path(OUT_DIR, case)
			dir.create(path = out_dir, showWarnings = FALSE)
			write.table(VariableFeatures(object), file = file.path(out_dir, "features.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
			object <- Dr(object, genes.list = VariableFeatures(object), dr = dr)
		}
	}
}

saveRDS(object, file = file.path(OUT_DIR, "object.Rds"))


sessionInfo()

q()


