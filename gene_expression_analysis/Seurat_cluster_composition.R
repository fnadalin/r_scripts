#!/usr/bin/R 

###### SEURAT V3 #######


ClustersComposition <- function(object, out.prefix, cl.ident.slot = "seurat_clusters", sample.ident.slot = "orig.ident") {

    sample_ident <- object@meta.data[[sample.ident.slot]]
    samples <- unique(sample_ident)

    cl_ident <- object@meta.data[[cl.ident.slot]]
    clusters <- unique(cl_ident)

    df <- matrix(NA, nrow = length(samples)*length(clusters), ncol = 3)
    df <- as.data.frame(df)
    colnames(df) <- c("sample", "cluster", "num")
    df$num <- 0 

    n <- 1
    for (s in rev(samples)) {
        for (c in clusters) {
            count <- sum(sample_ident == s & cl_ident == c)
            df[n,] <- c(s, c, count)
            n <- n + 1
        }
    }
    
    df$num <- as.numeric(df$num)
    df$cluster <- factor(df$cluster, levels = sort(clusters))
    write.table(df, file = paste0(out.prefix, ".tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

    g <- ggplot(df, aes(x=sample, y=num, fill=cluster)) + theme_minimal() + coord_flip() + xlab("") + theme(axis.text.y=element_text(size=15))
    g_abs <- g + geom_bar(stat="identity") + ylab("number of cells")
    pdf(paste0(out.prefix, "_SbyCabs.pdf"), width = 6, height = 3)
    print(g_abs)
    dev.off()

    g_norm <- g + geom_bar(stat="identity", position="fill") + ylab("fraction of cells")
    pdf(paste0(out.prefix, "_SbyCnorm.pdf"), width = 6, height = 3)
    print(g_norm)
    dev.off()

    g <- ggplot(df, aes(x=cluster, y=num, fill=sample)) + theme_minimal() + coord_flip() + xlab("") + theme(axis.text.y=element_text(size=15))
    g_abs <- g + geom_bar(stat="identity") + ylab("number of cells")
    pdf(paste0(out.prefix, "_CbySabs.pdf"), width = 6, height = 3)
    print(g_abs)
    dev.off()

    g_norm <- g + geom_bar(stat="identity", position="fill") + ylab("fraction of cells")
    pdf(paste0(out.prefix, "_CbySnorm.pdf"), width = 6, height = 3)
    print(g_norm)
    dev.off()

    return()
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: <obj.Robj> <out_dir> <params>\n")
	cat("\n<obj.Robj>   input Seurat object created at step 1\n")
	cat("<out_dir>    output directory containing the results of clusters composition\n")
	cat("<params>     file including the values for the parameters, separated by \"=\"\n\n")
	q()
}

library("Seurat")
library("cluster")
library("parallelDist")
library("factoextra")

OBJECT <- args[1]
OUT_DIR <- args[2]
PARAMS <- args[3]

params <- as.matrix(read.table(PARAMS, sep="="))
cores <- as.numeric(unlist(strsplit(params[params[,1]=="cores", 2], split=",")))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))

object <- readRDS(OBJECT)

cases <- c()
for (ymin in disps)
	cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
	cases <- c(cases, paste("vst_top", nfeat, sep=""))
drs <- paste("pca", cases, sep="_")

for (ll in drs) 
	dir.create(file.path(OUT_DIR, ll), showWarnings = FALSE, recursive = TRUE)

i <- 1
for (dr in drs) {
	print(paste("dr =", dr))
	for (kk in k) {
		for (r in res) {
			cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
			out_prefix <- file.path(OUT_DIR, dr, cl_ident)
			ClustersComposition(object = object, out.prefix = out_prefix, cl.ident.slot = cl_ident, sample.ident.slot = "sample.name") 
		}
	}
	i <- i + 1
}


sessionInfo()
q()

