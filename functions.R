
library("Matrix")
library("R.utils")
library("Seurat")
library("RColorBrewer")
library("ggplot2")
library("grid")
library("ggrepel")
library("scales") # to emulate ggplot default color palette
library("factoextra") # for clusters validation
library("cluster") # for clusters validation
library("pbapply")
# library("scater") # for MAST use in Seurat latest version N.B.: it is not available anymore on R version 3.6.3
library("limma") # for alias2SymbolTable
# library("AUCell") # for cells filtering
library("parallelDist")
library("reticulate") # to import python modules
library("future")
library("org.Hs.eg.db")
library("clusterProfiler")
library("ReactomePA")
library("plyr")
library("dplyr")
library("msigdbr")

MINIMUM <- 1e-300


ParseFeatureBarcodeMatrix <- function(matrix.dir, M) {

    mtx <- file.path(matrix.dir, "matrix.mtx")
    mtx_gz <- paste0(mtx, ".gz")

    if (!xor( file.exists(mtx), file.exists(mtx_gz) )) {
        write(paste("Either", mtx, "or", mtx_gz, "must exist"), stderr())
        q()
    } 

    compress <- !file.exists(mtx)
    features <- file.path(matrix.dir, "features.tsv")
    barcodes <- file.path(matrix.dir, "barcodes.tsv")

    if (!compress) {
        MM <- readMM(mtx)
        rownames(MM) <- as.matrix(read.table(features, sep="\t"))[,2]
        colnames(MM) <- drop(as.matrix(read.table(barcodes)))
    } else { 
        MM <- readMM(gzfile(paste0(mtx,".gz")))
        gz <- gzfile(paste0(features, ".gz"))
        rownames(MM) <- as.matrix(read.table(gz, sep="\t"))[,2]
        gz <- gzfile(paste0(barcodes, ".gz"))
        colnames(MM) <- drop(as.matrix(read.table(gz)))
    }

    eval.parent(substitute(M <- MM))
}


ParseFeatureBarcodeMatrixCombined <- function(matrix.dir, M_genes, M_gbc) {

    mtx <- file.path(matrix.dir, "matrix.mtx")
    mtx_gz <- paste0(mtx, ".gz")

    if (!xor( file.exists(mtx), file.exists(mtx_gz) )) {
        write(paste("Either", mtx, "or", mtx_gz, "must exist"), stderr())
        q()
    } 

    compress <- !file.exists(mtx)
    features <- file.path(matrix.dir, "features.tsv")
    barcodes <- file.path(matrix.dir, "barcodes.tsv")

    if (!compress) {
        MM <- readMM(mtx)
        FF <- as.matrix(read.table(features, sep="\t"))
        rownames(MM) <- FF[,2]
        feature_type <- FF[,3]
        colnames(MM) <- drop(as.matrix(read.table(barcodes)))
    } else { 
        MM <- readMM(gzfile(paste0(mtx,".gz")))
        gz <- gzfile(paste0(features, ".gz"))
        FF <- as.matrix(read.table(gz, sep="\t"))
        rownames(MM) <- FF[,2]
        feature_type <- FF[,3]
        gz <- gzfile(paste0(barcodes, ".gz"))
        colnames(MM) <- drop(as.matrix(read.table(gz)))
    }

    MM_genes <- MM[feature_type == "Gene Expression",]
    MM_gbc <- MM[feature_type == "CRISPR Guide Capture",]

    eval.parent(substitute(M_genes <- MM_genes))
    eval.parent(substitute(M_gbc <- MM_gbc))
}


ParseFeatureBarcodeMatrixExtract <- function(matrix.dir, M, assay) {

    mtx <- file.path(matrix.dir, "matrix.mtx")
    mtx_gz <- paste0(mtx, ".gz")

    if (!xor( file.exists(mtx), file.exists(mtx_gz) )) {
        write(paste("Either", mtx, "or", mtx_gz, "must exist"), stderr())
        q()
    }

    compress <- !file.exists(mtx)
    features <- file.path(matrix.dir, "features.tsv")
    barcodes <- file.path(matrix.dir, "barcodes.tsv")

    if (!compress) {
        MM <- readMM(mtx)
        FF <- as.matrix(read.table(features, sep="\t"))
        rownames(MM) <- FF[,2]
        feature_type <- FF[,3]
        colnames(MM) <- drop(as.matrix(read.table(barcodes)))
    } else {
        MM <- readMM(gzfile(paste0(mtx,".gz")))
        gz <- gzfile(paste0(features, ".gz"))
        FF <- as.matrix(read.table(gz, sep="\t"))
        rownames(MM) <- FF[,2]
        feature_type <- FF[,3]
        gz <- gzfile(paste0(barcodes, ".gz"))
        colnames(MM) <- drop(as.matrix(read.table(gz)))
    }

    MM_feat <- MM[feature_type == assay,]

    eval.parent(substitute(M <- MM_feat))
}


PCA <- function(object, out.dir, pcs.compute = 20, only_genes = NULL, suffix = "allgenes", reduction.name = "", do.jack = FALSE, ident.slot = "orig.ident", col = NULL) {

    if (!file.exists(out.dir))
        dir.create(file.path(getwd(), out.dir))

    if (reduction.name == "")
        reduction.name <- paste0("pca_", suffix)

    cat("Run PCA...")

    # PCA
    sink(file = file.path(out.dir, paste0(reduction.name, ".txt")))
    object <- RunPCA(object, features = only_genes, npcs = pcs.compute)
    sink()

    cat("done\n")

    cat("Generate PCA dotplots...")

    if (is.null(col)) {
        ident <- unique(object@meta.data[[ident.slot]])
        col <- hue_pal()(length(ident))
    }

    # plot the components
    plot.title <- paste(out.dir, "/", reduction.name, ".pdf", sep="")
    pdf(plot.title, width=21, height=30)
    ncol <- 4
    nrow <- 7
    title <- reduction.name
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
    grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
    # stdev explained by the first i components, for i=1,...,pcs.compute
    stdev <- 0
    n_plots <- floor(pcs.compute/2)
    for (i in 1:n_plots) {
        stdev1 <- stdev + object@reductions$pca@stdev[2*i-1]
        stdev2 <- stdev1 + object@reductions$pca@stdev[2*i]
        title <- paste("stdev(pc", 2*i-1, ") = ", round(stdev1, digits=1), ", stdev(pc", 2*i, ")=", round(stdev2, digits=1), "\n", sep="")
        g <- DimPlot(object, dims = c(2*i-1, 2*i), cols = col, group.by = ident.slot) + ggtitle(title)
        print(g, vp = viewport(layout.pos.row = ((i-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (i-1)%%ncol+1))
        if (i%%(nrow*ncol) == 0 && i < n_plots) {
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
        }
        stdev <- stdev2
    }
    dev.off()

    cat("done\n")

    ### NB: The functions used below have reduction.type hardcoded ("pca")
    # So the correct pca ID is assigned at the end

    if (do.jack) {
        cat("Generate PCA JackStraw plot...")
    
        # Significant components
        object <- JackStraw(object, dims = pcs.compute, num.replicate = 100)
        plot.title <- paste(out.dir, "/Jack_", reduction.name, ".pdf", sep="")
        pdf(plot.title, width = 7, height = 7*ceiling(pcs.compute/25))
        print(JackStrawPlot(object, nCol = 5, PCs = 1:pcs.compute))
        dev.off()

        cat("done\n")
    }

    cat("Generate PCA elbow plot...")

    # PCA stdev plot
    pdf(paste(out.dir, "/stdev_", reduction.name, ".pdf", sep=""))
    print(ElbowPlot(object, ndims = pcs.compute))
    dev.off()

    cat("done\n")

    names(object@dr)[names(object@dr) == "pca"] <- reduction.name

    return(object)
}


CellsClusters <- function(object, dr = "pca", dim = 10, k = c(30), res = 0.1*(1:10)) {

    cat("Cluster cells from PCA...\n")

    if (!any(names(object@reductions) == dr)) {
        cat(paste("Reduction ", dr, " has not been computed yet\n", sep=""))
        return(object)
    }

    for (kk in k) {
        object <- FindNeighbors(object, k.param = kk, force.recalc = TRUE, reduction = dr, dims = 1:dim)
        for (r in res) {
            cl_ident <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            object <- FindClusters(object, resolution = r)
            names(object@meta.data)[names(object@meta.data) == "seurat_clusters"] <- cl_ident
            object@active.ident <- as.factor(object$orig.ident)
        }
    }

    object <- SetIdent(object, value = "orig.ident")
    return(object)
}


CellsClustersLeiden <- function(object, dr = "pca", dim = 10, k = c(30), res = 0.1*(1:10)) {

    cat("Cluster cells from PCA...\n")

    if (!any(names(object@reductions) == dr)) {
        cat(paste("Reduction ", dr, " has not been computed yet\n", sep=""))
        return(object)
    }

    for (kk in k) {
        object <- FindNeighbors(object, k.param = kk, force.recalc = TRUE, reduction = dr, dims = 1:dim)
        for (r in res) {
            cl_ident <- paste0("clusters_leiden_", dr, dim, "_k", kk, "_res", r)
            object <- FindClusters(object, resolution = r, algorithm = 4, method = "igraph")
            names(object@meta.data)[names(object@meta.data) == "seurat_clusters"] <- cl_ident
            object@active.ident <- as.factor(object$orig.ident)
        }
    }

    object <- SetIdent(object, value = "orig.ident")
    return(object)
}


Silhouette <- function(object, out.dir, dr = "pca", dim = 10, k = 30, res = 0.1*(1:10)) {

    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    data <- as.matrix(object@reductions[[dr]]@cell.embeddings[,1:dim])
    dist <- parDist(data, threads = THREADS)
    for (kk in k) {
        for (r in res) {
            cl_ident <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            out_prefix <- file.path(out.dir, cl_ident)
            ident = c(as.matrix(object@meta.data[cl_ident]))
            id <- unique(ident)
            if (length(id) > 1) { # FIXED: without this check, it fails if there is only one cluster
                x <- plyr::mapvalues(ident, from = id, to = 1:length(id))
                silh <- silhouette(as.numeric(x), dist)
                sink(paste(out_prefix, ".txt", sep=""))
                print(summary(silh))
                sink()
                pdf(paste(out_prefix, ".pdf", sep=""))
                plot(silh)
                dev.off()
            }
        }
    }

    return()
}


SilhouetteNEW <- function(object, out.dir, dr = "pca", dim = 10, k = c(30), res = 0.1*(1:10)) {

    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    data <- as.matrix(object@reductions[[dr]]@cell.embeddings[,1:dim])
    dist <- parDist(data, threads = THREADS)
    for (kk in k) {
        for (r in res) {
            cl_ident <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            out_prefix <- file.path(out.dir, cl_ident)
            ident = c(as.matrix(object@meta.data[cl_ident]))
            id <- sort(unique(ident))
            if (length(id) > 1) { # FIXED: without this check, it fails if there is only one cluster
                x <- plyr::mapvalues(ident, from = id, to = 1:length(id))
                silh <- silhouette(as.numeric(x), dist)
                sink(paste(out_prefix, ".txt", sep=""))
                print(summary(silh))
                sink()
                pdf(paste(out_prefix, ".pdf", sep=""))
                print(fviz_silhouette(silh))
                dev.off()
            }
        }
    }

    return()
}


SamplesClustersComposition <- function(object, out.dir, dr = "pca", dim = 10, k = 30, res = 0.1*(1:10), sample.ident.slot = "orig.ident") {

    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    sample_ident <- object@meta.data[[sample.ident.slot]]
    samples <- unique(sample_ident)

    for (kk in k) {
        for (r in res) {
            cl.ident.slot <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            cl_ident <- object@meta.data[[cl.ident.slot]]
            clusters <- unique(cl_ident)
            out_prefix <- file.path(out.dir, cl_ident)

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
            write.table(df, file = file.path(out.dir, paste0(cl.ident.slot, ".tsv")), row.names=FALSE, sep = "\t", quote=FALSE)

            g <- ggplot(df, aes(x=sample, y=num, fill=cluster)) + theme_minimal() + coord_flip() + xlab("") + theme(axis.text.y=element_text(size=15))
            g_abs <- g + geom_bar(stat="identity") + ylab("number of cells")
            pdf(file.path(out.dir, paste0(cl.ident.slot, "_SbyCabs.pdf")), width=6, height=3)
            print(g_abs)
            dev.off()

            g_norm <- g + geom_bar(stat="identity", position="fill") + ylab("fraction of cells")
            pdf(file.path(out.dir, paste0(cl.ident.slot, "_SbyCnorm.pdf")), width=6, height=3)
            print(g_norm)
            dev.off()

            g <- ggplot(df, aes(x=cluster, y=num, fill=sample)) + theme_minimal() + coord_flip() + xlab("") + theme(axis.text.y=element_text(size=15))
            g_abs <- g + geom_bar(stat="identity") + ylab("number of cells")
            pdf(file.path(out.dir, paste0(cl.ident.slot, "_CbySabs.pdf")), width=6, height=3)
            print(g_abs)
            dev.off()

            g_norm <- g + geom_bar(stat="identity", position="fill") + ylab("fraction of cells")
            pdf(file.path(out.dir, paste0(cl.ident.slot, "_CbySnorm.pdf")), width=6, height=3)
            print(g_norm)
            dev.off()
        }
    }

    return()
}


PlotTSNEClusters <- function(object, out.file, dr = "pca", dim = 10, k = c(30), res = 0.1*(1:10)) {

    pdf(out.file, width = 8.27, height = 11.69)
    ncol <- 2
    nrow <- 3
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol, heights = unit(rep(1/nrow, nrow), "npc"))))
    i <- 1
    for (kk in k) {
        for (r in res) {
            title <- cl.ident.slot <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            g <- DimPlot(object, reduction = "tsne", group.by = cl.ident.slot, label = FALSE) + ggtitle(title)
            print(g, vp = viewport(layout.pos.row = ((i-1)%%(ncol*nrow))%/%ncol+1, layout.pos.col = (i-1)%%ncol+1))
            if (i%%(nrow*ncol) == 0) {
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
            }
            i <- i + 1
        }
    }
    dev.off()

    return()
}


PlotPCAClusters <- function(object, out.file, dr = "pca", dim = 10, k = c(30), res = 0.1*(1:10)) {

    pdf(out.file, width = 8.27, height = 11.69)
    ncol <- 2
    nrow <- 3
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol, heights = unit(rep(1/nrow, nrow), "npc"))))
    i <- 1
    for (kk in k) {
        for (r in res) {
            cl.ident.slot <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            title <- paste0("clusters_", dr, dim, "_k", kk, "_res", r)
            g <- DimPlot(object, reduction = "pca", group.by = cl.ident.slot, label = FALSE) + ggtitle(title)
            print(g, vp = viewport(layout.pos.row = ((i-1)%%(ncol*nrow))%/%ncol+1, layout.pos.col = (i-1)%%ncol+1))
            if (i%%(nrow*ncol) == 0) {
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
            }
            i <- i + 1
        }
    }
    dev.off()

    return()
}


# the logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
# specify the assay explicitely in order to make sure that the integrated data will NOT be used for differential expression
GeneMarkersTable <- function(object, out.name, ident.1, ident.2, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

    cells.1 <- WhichCells(object = object, idents = ident.1)
    exp1 <- apply(GetAssayData(object = object, assay = "RNA")[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

    cells.2 <- WhichCells(object = object, idents = ident.2)
    exp2 <- apply(GetAssayData(object = object, assay = "RNA")[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

    cluster_markers <- FindMarkers(object = object, assay = "RNA", ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use)
    df <- data.frame(geneID = rownames(cluster_markers),
                     pct.1 = cluster_markers$pct.1, pct.2 = cluster_markers$pct.2,
                     avg_exp.1 = exp1[rownames(cluster_markers)], avg_exp.2 = exp2[rownames(cluster_markers)],
                     avg_log2FC = cluster_markers$avg_log2FC,
                     p_val = cluster_markers$p_val,
                     p_val_adj = cluster_markers$p_val_adj)
        
    write.table(df, file=out.name, sep="\t", quote=FALSE, row.names = FALSE)

    return()
}


# the logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
GeneMarkersTableNEW <- function(object, out.name, ident.1, ident.2, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

    cells.1 <- WhichCells(object = object, idents = ident.1)
    exp1 <- apply(GetAssayData(object = object, assay = "RNA")[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

    cells.2 <- WhichCells(object = object, idents = ident.2)
    exp2 <- apply(GetAssayData(object = object, assay = "RNA")[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

    # NEW: specify the assay explicitely!!!
    # this is to make sure that the integrated data will NOT be used for differential expression
    cluster_markers <- FindMarkers(object = object, assay = "RNA", ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use)
    df <- data.frame(geneID = rownames(cluster_markers),
                     pct.1 = cluster_markers$pct.1, pct.2 = cluster_markers$pct.2,
                     avg_exp.1 = exp1[rownames(cluster_markers)], avg_exp.2 = exp2[rownames(cluster_markers)],
                     avg_log2FC = cluster_markers$avg_log2FC,
                     p_val = cluster_markers$p_val,
                     p_val_adj = cluster_markers$p_val_adj)
    
    write.table(df, file=out.name, sep="\t", quote=FALSE, row.names = FALSE)

    return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
# NEW TABLE FORMAT
ClusterGeneMarkersByPairs <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

    if (id != "") Idents(object) <- id
    if (length(clusters) == 0) clusters <- levels(object@active.ident)

    for (i in 1:(length(clusters)-1)) {    
        for (j in (i+1):length(clusters)) {
            deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv", sep="")
            GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = clusters[j], test.use = test.use, min.pct = min.pct, logFC = logFC)            
        }
    }

    return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
# NEW TABLE FORMAT
ClusterGeneMarkersVsAll <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

    if (id != "") Idents(object) <- id
    if (length(clusters) == 0) clusters <- levels(object@active.ident)

    for (i in 1:length(clusters)) {
        deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
        GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = setdiff(clusters, c(clusters[i])), test.use = test.use, min.pct = min.pct, logFC = logFC)    
    }

    return()
}


FilterDEGtable <- function(table, logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, sort = FALSE, genes.use = NULL, genes.filter = NULL) {

    # filter out genes based on abs(logFC), adj. pval, and percentage of cells where the gene is detected
    # IMPORTANT!!! Transform to base 2 log
    v <- table[which(abs(table$avg_log2FC) > logFC.filt & table$p_val_adj < adjpval.filt & (table$pct.1 > min.pct | table$pct.2 > min.pct)),]

    # exclude undefined p-values and HIV genes
    v <- v[!is.na(v[,1]),]

    # exclude gene.filter genes / keep only specified genes
    if (!is.null(genes.filter)) v <- v[!(v$geneID %in% genes.filter),]
    if (!is.null(genes.use)) v <- v[v$geneID %in% genes.use,]

    # print the top num DEG, according to the value of ad. p-value
    num <- min(num, nrow(v))
    v <- v[1:num,]
    if (sort) {
        v <- v[order(v$avg_log2FC, decreasing=TRUE),]
    }
    
    return(v)
}


FilterClusterGeneMarkers <- function(input.table, output.table, logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, sort = FALSE, genes.use = NULL, genes.filter = NULL) {

    m <- read.table(input.table, sep="\t", header=TRUE)
    v <- FilterDEGtable(table = m, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, sort = sort, genes.use = genes.use, genes.filter = NULL)
    write.table(v, file=output.table,  sep="\t", quote=FALSE, row.names = FALSE)
    
    return()
}


HeatmapPlot <- function(object, filt.table, plot.name, cells = NULL, width = 5, height = 7) {
        
    table <- read.table(filt.table, sep="\t", header = TRUE)
    genes <- as.character(table$geneID)
    if (length(genes) == 0) next

    pdf(plot.name, width = width, height = height)    
    print(DoHeatmap(object = object, cells = cells, features = genes, label = TRUE, size = 2))
    dev.off()

    return()
}


HeatmapPlotNEW <- function(object, filt.table, plot.name, cells = NULL, width = 5, height = 7) {
        
    table <- read.table(filt.table, sep="\t", header = TRUE)
    genes <- as.character(table$geneID)
    if (length(genes) == 0) next

    pdf(plot.name, width = width, height = height)    
    print(DoHeatmap(object = object, cells = cells, features = genes, label = TRUE, size = 2, assay = "RNA"))
    dev.off()

    return()
}


# Filt table can just be a list of gene names
# color.genes: array of genes to be plot in color
# label.genes: array of genes to be labelled
# VolcanoPlotNEW <- function(all.table, filt.table, plot.name, title = "", revert = FALSE, genes.use = NULL) {
VolcanoPlot <- function(all.table, color.genes = NULL, label.genes = NULL, plot.name, title = "", subtitle = "", group.name = "top50sign",
                        revert = FALSE, genes.use = NULL, genes.filter = NULL, colors = c("gray", "blue"), point.size = 1) {

    all <- read.table(all.table, sep="\t", header=TRUE)
    if (!is.null(genes.use)) all <- all[all$geneID %in% genes.use,]
    if (!is.null(genes.filter)) all <- all[!(all$geneID %in% genes.filter),]

    if (revert) all$avg_log2FC <- -all$avg_log2FC

    maxFC <- max(all$avg_log2FC[which(is.finite(all$avg_log2FC))])
    minFC <- min(all$avg_log2FC[which(is.finite(all$avg_log2FC))])
    
    non_detectable <- which(all$p_val_adj < MINIMUM)
    if (length(non_detectable) > 0) 
        all$p_val_adj[non_detectable] <- rep(MINIMUM, length(non_detectable))

    inf_pos <- which(!is.finite(all$avg_log2FC) && all$avg_log2FC > 0)
    inf_neg <- which(!is.finite(all$avg_log2FC) && all$avg_log2FC < 0)
    if (length(inf_pos) > 0) 
        all$avg_log2FC[inf_pos] <- rep(minFC,length(inf_pos))
    if (length(inf_neg) > 0) 
        all$avg_log2FC[inf_neg] <- rep(minFC,length(inf_neg))

    df <- data.frame(x = all$avg_log2FC, y = -log10(all$p_val_adj), z = all$geneID)
    df$type <- factor(ifelse(all$geneID %in% color.genes, group.name, "nogroup"))
    df$type <- factor(df$type, levels = c("nogroup", group.name))
    df$col <- factor(ifelse(all$geneID %in% color.genes, "col", "notcol"))
    df$lab <- factor(ifelse(all$geneID %in% label.genes, "lab", "notlab"))

    pdf(plot.name, width=6.5, height=7)

    maxabsFC <- max(abs(minFC), abs(maxFC))

    g <- ggplot(data = df, aes(x = x, y = y, color=type)) + theme_bw() + xlab(expression(log[2](FC))) + ylab(expression(-log[10](FDR))) + xlim(-maxabsFC, maxabsFC)
    g <- g + geom_point(size=point.size)
    g <- g + scale_color_manual(values=colors)
    g <- g + geom_text_repel(data = df[which(df$lab == "lab"),], aes(label = z), color="black")
    g <- g + ggtitle(label=title, subtitle = subtitle) + theme(plot.title=element_text(size=20), plot.subtitle=element_text(size=10))
    g <- g + theme(text=element_text(size=15)) + theme(legend.position = "none")
    print(g)
    
    dev.off()

    return()
}


VolcanoPlotFilter <- function(all.table, filt.table, plot.name, title = "", subtitle = "", genes.use = NULL, genes.filter = NULL, revert = FALSE) {

    filt <- read.table(filt.table, sep="\t", header=TRUE)

    label.genes <- color.genes <- as.character(filt$geneID) # otherwise they are seen as factors!!

    VolcanoPlot(all.table = all.table, color.genes = color.genes, label.genes = label.genes, 
        plot.name = plot.name, title = title, subtitle = subtitle, genes.use = genes.use, genes.filter = genes.filter, revert = revert)

    return()
}


FilterClusterGeneMarkersPairs <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, genes.use = NULL, heatmap = FALSE, volcano = FALSE) {

    if (id != "") Idents(object) <- id
    if (length(clusters) == 0) clusters <- levels(Idents(object))

    dir.create(paste(out.dir, "heatmap", sep="/"), showWarnings = FALSE)
    dir.create(paste(out.dir, "volcano", sep="/"), showWarnings = FALSE)

    for (i in 1:(length(clusters)-1)) {
        for (j in (i+1):length(clusters)) {
            deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv", sep="")
            filt.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "-filtered.tsv", sep="")
            if (file.exists(deg.table)) {
                FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, genes.use = genes.use)
                if (heatmap) {
                    cells <- colnames(GetAssayData(object))[Idents(object) %in% sort(c(clusters[i], clusters[j]))]
                    plot.name <- paste(out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".pdf", sep="")
                    HeatmapPlot(object, filt.table, plot.name, cells = cells)
                }
                if (volcano) {
                    volcano.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".pdf", sep="")
                    title <- paste("cl", clusters[i], "-", clusters[j], sep="")
                    VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano.plot, title = title, genes.use = genes.use)

                    volcano_rev.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[j], "-", clusters[i], ".pdf", sep="")
                    title_rev <- paste("cl", clusters[j], "-", clusters[i], sep="")
                    VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano_rev.plot, title = title_rev, revert = TRUE, genes.use = genes.use)
                }
            }
        }
    }
    
    return()
}

FilterClusterGeneMarkersAll <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, genes.use = NULL, heatmap = FALSE, volcano = FALSE) {

    dir.create(paste(out.dir, "heatmap", sep="/"))
    dir.create(paste(out.dir, "volcano", sep="/"))

    if (id != "") Idents(object) <- id
    if (length(clusters) == 0) clusters <- levels(Idents(object))
    cells <- colnames(GetAssayData(object))[Idents(object) %in% clusters]

    for (i in 1:length(clusters)) {
        deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
        filt.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all-filtered.tsv", sep="")
        if (file.exists(deg.table)) {
        FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, genes.use = genes.use)
            if (heatmap) {
                plot.name <- paste(out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-all.pdf", sep="")
                HeatmapPlot(object, filt.table, plot.name, cells = cells, width = 9, height = 7)
            }
            if (volcano) {
                volcano.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-all.pdf", sep="")
                title <- paste("cl", clusters[i], "-all", sep="")
                VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano.plot, title = title, genes.use = genes.use)
            }
        }
    }
    
    return()
}


# remove the cells that belong to the specified clusters
# clusters are specified as a comma-separated string
# return the list of selected cells 
# select = TRUE will select the specified identities instead of removing them
FilterOutClusters <- function (object, id, cl_filter, select = FALSE) {

    if (cl_filter == "") 
        return(colnames(GetAssayData(object)))

    cl <- unlist(strsplit(cl_filter, split=","))    

    if (select) {
        filtered_cells <- names(object@active.ident[object@meta.data[[id]] %in% cl])
    } else {
        filtered_cells <- names(object@active.ident[!(object@meta.data[[id]] %in% cl)])
    }

    return(filtered_cells)
}


# reassign cluster names
# renaming is specified as a comma-separated string and refers to clusters 0,1,2... in this order
# specify the array of clusters to be masked (i.e. exlcuded from the renaming
# TODO: add a new column to the metadata that contains the renamed clusters!!!
RenameClusters <- function (object, renaming, id, cl_filter = NULL) {

    old.cl.names <- orig.cl.names <- as.character(sort(as.numeric(unique(object@meta.data[[id]]))))
    if (!is.null(cl_filter)) {
        cl <- unlist(strsplit(cl_filter, split=","))
        old.cl.names <- orig.cl.names[!(orig.cl.names %in% cl)]
        object@meta.data[[id]] <- plyr::mapvalues(x = object@meta.data[[id]], from = cl, to = rep("X", length(cl))) # assign X to removed clusters
    }
    new.cl.names <- unlist(strsplit(renaming, split=","))
    object@meta.data[[id]] <- plyr::mapvalues(x = object@meta.data[[id]], from = old.cl.names, to = new.cl.names) # rename only filtered clusters according to the mapping

    return(object)
}


SelectUpDown <- function (table, revert = FALSE) {

    if (!revert) {
        up <- "up"
        down <- "down"
    } else {
        up <- "down"
        down <- "up"
    }

    g <- rep(up, nrow(table))
    g[which(table$avg_log2FC < 0)] <- down
    df <- data.frame(geneID = table$geneID,    group = g)
    
    return(df)
}

GOenrichPreprocessing <- function(deg.table, sel.degs, org, min.perc = 0.1, logFC.filt = 1, adjpval.filt = 0.1, revert = FALSE, num = 1000000, to.official = FALSE) {

    m <- read.table(deg.table, sep="\t", header=TRUE)
    t <- FilterDEGtable(table = m, min.pct = min.perc, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, num = num)
    df <- SelectUpDown(table = t, revert = revert)
    
    orglib <- paste0("org.", org, ".eg.db")

    # convert to official gene symbol
    geneID <- geneIDorig <- as.character(df$geneID)
    if (to.official) 
        geneID <- alias2SymbolTable(geneIDorig, species = org) # to convert gene aliases to official gene names

    # convert to Entrex ID
    # NB: redundant terms are collapsed, this includes NA as well
    SymbToEnt <- bitr(geneID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orglib, drop = FALSE)

    # create the df with the 3 gene nomenclatures
    df <- data.frame(Entrez = rep(NA,length(geneID)), geneID = geneID, geneIDorig = geneIDorig, group = df$group)
    df$Entrez[match(SymbToEnt$SYMBOL, df$geneID)] <- SymbToEnt$ENTREZID
    write.table(df, file = sel.degs, sep="\t", quote = FALSE, row.names = FALSE)

    return()
}


# species is mandatory
# keep track of the original gene symbol (which may be different from the geneID converted to entrez if to.official == TRUE)
# print the data frame containing the gene symbol, the original gene symbol, the entrez ID, and the group (up or down)
GOenrichNEW <- function(table, out_prefix, org, types = c("MF", "BP", "CC"), pval = 1, qval = 0.05, level = 4, similarity = 0.7, cat.to.show = 15, title = "") {

    orglib <- paste("org.", org, ".eg.db", sep="")

    df <- read.table(table, header = TRUE)
    df$group <- factor(df$group, levels = c("up", "down"))

    # create the input df for pw enrichment analysis
    # it only contains unique non-NA Entrez IDs
    mapped <- which(!is.na(df$Entrez))
    df_map <- df[mapped,]
    df_map_nr <- unique(data.frame(Entrez = df_map$Entrez, group = df_map$group))

    for (type in types) {

        # perform GO enrichment analysis
        ego <- compareCluster(Entrez~group, data=df_map_nr, fun="enrichGO", OrgDb=orglib, ont=type, pAdjustMethod="BH", pvalueCutoff=pval, qvalueCutoff=qval, readable=TRUE)

        if (length(ego@compareClusterResult$Description)  > 0) {
    
            # set the minimum granularity level
            # ego <- gofilter(ego, level = LEVEL)
            # remove redundant GO terms and by keeping a representative terms
            # the representative term is the one with lowest p-value
            ego <- simplify(ego, cutoff = similarity, by = "p.adjust", select_fun = min)

            write.table(as.data.frame(ego), file = paste0(out_prefix, "_", type, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE)
        
            PLOT <- paste0(out_prefix, "_", type, ".pdf")
            height <- max(4, 0.5+0.4*min(length(ego@compareClusterResult$Description), cat.to.show))
            plot.title <- ""
            if (title != "") plot.title <- paste0(title, " - ")
            plot.title <- paste0(plot.title, "GO enrichment (", type, ")")
            pdf(PLOT, width = 12, height = height)
            print(dotplot(ego, showCategory = cat.to.show, title = plot.title))
            dev.off()
        }
    }

    return()
}


# db is mandatory: "reactome" or "kegg"
# species is mandatory: "Hs", "Mm", ...
# orgname is mandatory: "human", "mouse", ...
# keep track of the original gene symbol (which may be different from the geneID converted to entrez if to.official == TRUE)
# print the data frame containing the gene symbol, the original gene symbol, the entrez ID, and the group (up or down)
PWenrichNEW <- function(table, out_prefix, db, org, orgname, pval = 1, qval = 0.05, cat.to.show = 15, title = "") {

    orglib <- paste0("org.", org, ".eg.db")

    df <- read.table(table, header = TRUE)
    df$group <- factor(df$group, levels = c("up", "down"))

    # create the input df for pw enrichment analysis
    # it only contains unique non-NA Entrez IDs
    mapped <- which(!is.na(df$Entrez))
    df_map <- df[mapped,]
    df_map_nr <- unique(data.frame(Entrez = df_map$Entrez, group = df_map$group))

    if (db == "reactome") 
        pwe <- compareCluster(Entrez~group, data=df_map_nr, fun="enrichPathway", organism=orgname, pvalueCutoff=pval, pAdjustMethod="BH", qvalueCutoff=qval, readable=TRUE)
    if (db == "kegg") 
        pwe <- compareCluster(Entrez~group, data=df_map_nr, fun="enrichKEGG", organism=orgname, pvalueCutoff=pval, pAdjustMethod="BH", qvalueCutoff=qval)
    
    if (length(pwe@compareClusterResult$Description) > 0) {

        df <- as.data.frame(pwe)

        # compute the sample odds ratio
        odds_ratio <- apply(df[,c("GeneRatio","BgRatio")], 1, function(x) {
            set <- as.numeric(unlist(strsplit(x[1], split = "/")))
            comp <- as.numeric(unlist(strsplit(x[2], split = "/")))
            (set[1]/(comp[1]-set[1])) / ((set[2]-set[1])/(comp[2]-set[2]-comp[1]+set[1])) }
        )

        df_with_odds <- cbind(df, odds_ratio = odds_ratio)
        write.table(df_with_odds, file = paste0(out_prefix, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)

        PLOT <- paste0(out_prefix, ".pdf")
        height <- max(4,0.5+0.4*min(length(pwe@compareClusterResult$Description), cat.to.show))
        if (title != "") title <- paste(title, " - ", sep=" ")
        title <- paste0(title, "PW enrichment (", db, ")")
        pdf(PLOT, width = 12, height = height)
        print(dotplot(pwe, showCategory = cat.to.show, title = title))
        dev.off()
    }

    return()
}


# org is mandatory
KEGGpostprocessingNEW <- function(input.table, output.table, org) {

    data <- read.table(input.table, header=TRUE, sep="\t")

    data_conv <- data
    idx <- which(colnames(data_conv) == "geneID")
    data_conv$geneID <- apply(data, 1, function(x) {
        v <- unlist(strsplit(as.character(x[idx]), split="/"))
        Ent2Symb <- bitr(v, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = paste("org", org, "eg.db", sep="."), drop=FALSE)
        return(paste(Ent2Symb$SYMBOL, collapse="/"))
        }
    )

    write.table(data_conv, file=output.table, row.names=FALSE, sep="\t", quote=FALSE)

    return()
}


# species is mandatory
# keep track of the original gene symbol (which may be different from the geneID converted to entrez if to.official == TRUE)
# print the data frame containing the gene symbol, the original gene symbol, the entrez ID, and the group (up or down)
MSDenrich <- function(table, out_prefix, species, org, types = c("H"), pval = 1, qval = 0.05, cat.to.show = 15, title = "") {

    orglib <- paste("org.", org, ".eg.db", sep="")

    df <- read.table(table, header = TRUE)
    df$group <- factor(df$group, levels = c("up", "down"))

    # create the input df for pw enrichment analysis
    # it only contains unique non-NA Entrez IDs
    mapped <- which(!is.na(df$Entrez))
    df_map <- df[mapped,]
    df_map_nr <- unique(data.frame(Entrez = df_map$Entrez, group = df_map$group))

    for (type in types) {

        m_t2g <- msigdbr(species = species, category = type) %>% dplyr::select(gs_name, entrez_gene)

        # perform GO enrichment analysis
        ego <- compareCluster(Entrez~group, data=df_map_nr, fun="enricher", TERM2GENE=m_t2g, pAdjustMethod="BH", pvalueCutoff=pval, qvalueCutoff=qval)

        if (length(ego@compareClusterResult$Description)  > 0) {

            write.table(as.data.frame(ego), file = paste0(out_prefix, "_", type, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE)
        
            PLOT <- paste0(out_prefix, "_", type, ".pdf")
            height <- max(4, 0.5+0.4*min(length(ego@compareClusterResult$Description), cat.to.show))
            plot.title <- ""
            if (title != "") plot.title <- paste0(title, " - ")
            plot.title <- paste0(plot.title, "MSD enrichment (", type, ")")
            pdf(PLOT, width = 12, height = height)
            print(dotplot(ego, showCategory = cat.to.show, title = plot.title))
            dev.off()
        }
    }

    return()
}




