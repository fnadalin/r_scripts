
# Input:
#  - list of directories with the results for all samples
#  - list of sample names
#  - list of sample names for the samples containing the surviving clones
# Output:
#  - table with the union of the surviving clones and the list of samples where they are found


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--results_dir_list", type = "character",
        help="[REQUIRED] tsv file containing the files with the CB classification for each sample, one per line, where each column represents:\n		1 - the sample name\n		2 - the directory containing the results of the CB classification\n		3 - whether the sample contains surviving clones (\"surv\") or not (\"non-surv\")"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] output dir where the results are stored")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$results_dir_list)) {
    write("Option --results_dir_list is required\nTry --help for help", stderr()) 
    q()
} else {
    RESULTS_DIR_LIST <- opt$results_dir_list
}

if (is.null(opt$out_dir)) {
    write("Option --out_dir is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_DIR <- opt$out_dir
}


################################### IMPORT #####################################

library("funr")

SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())

# print parameters to file
sink(file = file.path(OUT_DIR, paste0("parameters_", SCRIPT_NAME)))
for (i in 1:length(opt))
    cat(paste0(names(opt)[i], "=", opt[[i]], "\n"))
sink()


################################## EXECUTION ###################################

# output file names

library("VennDiagram")
library("ggplot2")

if (!dir.exists(OUT_DIR))
    dir.create(OUT_DIR, recursive = TRUE)

IN_CLONE_INFO              <- "stats_clones.tsv"
IN_CELL_INFO               <- "CB_classification.tsv"

OUT_CELL_INFO              <- "CB_classification_surv.tsv"
OUT_TABLE_CLONE_SAMPLE     <- file.path(OUT_DIR, "clone_sample_CB_count.tsv")
OUT_TABLE_CLONE_SAMPLE_SINGLE <- file.path(OUT_DIR, "clone_sample_CB_count_single.tsv")
OUT_TABLE_CLONE_SAMPLE_SINGLE_AND_DOUBLET <- file.path(OUT_DIR, "clone_sample_CB_count_single_and_doublet.tsv")
OUT_TABLE_CLONE_SAMPLE_COINFECTED <- file.path(OUT_DIR, "clone_sample_CB_count_coinfected.tsv")
OUT_TABLE_CELL_CLASS       <- file.path(OUT_DIR, "barplot_CB_classification.tsv")
OUT_TABLE_CLONE_NUM        <- file.path(OUT_DIR, "barplot_GBC_count.tsv")
OUT_TABLE_CELL_SURV_CLASS  <- file.path(OUT_DIR, "barplot_CB_surv_classification.tsv")
OUT_TABLE_SURV_CLONES_UNION <- file.path(OUT_DIR, "surv_clones_union.tsv")
OUT_TABLE_SURV_CLONES_INTERSECTION <- file.path(OUT_DIR, "surv_clones_intersection.tsv")
OUT_FIGURE_CELL_CLASS      <- file.path(OUT_DIR, "barplot_CB_classification.pdf")
OUT_FIGURE_CELL_CLASS_FRAC <- file.path(OUT_DIR, "barplot_CB_classification_frac.pdf")
OUT_FIGURE_CLONE_NUM       <- file.path(OUT_DIR, "barplot_GBC_count.pdf")
OUT_FIGURE_CELL_SURV_CLASS <- file.path(OUT_DIR, "barplot_CB_surv_classification.pdf")
OUT_FIGURE_CLONE_VENN      <- file.path(OUT_DIR, "clone_sample_early_venn.png")
OUT_FIGURE_CLONE_SURV_VENN <- file.path(OUT_DIR, "clone_sample_surv_venn.png")
OUT_FIGURE_CLONE_SINGLE_VENN      <- file.path(OUT_DIR, "clone_sample_early_single_venn.png")
OUT_FIGURE_CLONE_SURV_SINGLE_VENN <- file.path(OUT_DIR, "clone_sample_surv_single_venn.png")
OUT_FIGURE_CLONE_SINGLE_AND_DOUBLET_VENN <- file.path(OUT_DIR, "clone_sample_early_single_and_doublet_venn.png")
OUT_FIGURE_CLONE_SURV_SINGLE_AND_DOUBLET_VENN <- file.path(OUT_DIR, "clone_sample_surv_single_and_doublet_venn.png")
OUT_FIGURE_CLONE_COINF_VENN      <- file.path(OUT_DIR, "clone_sample_early_coinf_venn.png")
OUT_FIGURE_CLONE_SURV_COINF_VENN <- file.path(OUT_DIR, "clone_sample_surv_coinf_venn.png")
OUT_FIGURE_CLONE_SURV_CORR <- file.path(OUT_DIR, "clone_sample_surv_correlation.pdf")
OUT_FIGURE_HIST_CLONE_BY_CELL_NUM <- file.path(OUT_DIR, "hist_early_samples_clones_by_cell_num.pdf")

# parse the input

results_dir_list <- as.matrix(read.table(RESULTS_DIR_LIST, sep = "\t"))
sample_name <- results_dir_list[,1]
sample_dir <- results_dir_list[,2]
sample_surv_name <- sample_name[results_dir_list[,3] == "surv"]

idx_surv <- which(sample_name %in% sample_surv_name)
idx_before_surv <- which(!(sample_name %in% sample_surv_name))

n_samples <- nrow(results_dir_list)

# collect CB classification results

cat("collect CB classification results\n")

# generate the barplot with the number of cells, grouped by class
df_cell_bar <- as.data.frame(matrix(NA, nrow = 0, ncol = 3))
colnames(df_cell_bar) <- c("sample", "class", "count")
for (i in 1:n_samples) {
    filename <- file.path(sample_dir[i], IN_CELL_INFO)
    df <- read.table(filename, sep = "\t")
    t <- table(df$class)
    df <- data.frame(sample = rep(sample_name[i], length(t)), t)
    colnames(df) <- colnames(df_cell_bar)
    df_cell_bar <- rbind(df_cell_bar, df)
}
df_cell_bar$class <- factor(df_cell_bar$class, levels = rev(c("single_guide", "coinfected", "doublet", "uninfected")))

write.table(df_cell_bar, file = OUT_TABLE_CELL_CLASS, sep = "\t", row.names = FALSE, quote = FALSE)

g <- ggplot(data = df_cell_bar, aes(x = sample, y = count, fill = class)) + theme_classic() + geom_bar(stat = "identity") 
g <- g + xlab("") + ylab("Number of CB") + ggtitle("CB classification")
pdf(OUT_FIGURE_CELL_CLASS, width = 4.5, height = 3)
g
dev.off()

g <- ggplot(data = df_cell_bar, aes(x = sample, y = count, fill = class)) + theme_classic() + geom_bar(stat = "identity", position = "fill") 
g <- g + xlab("") + ylab("Number of CB") + ggtitle("CB classification")
pdf(OUT_FIGURE_CELL_CLASS_FRAC, width = 4.5, height = 3)
g
dev.off()

# create a clone / sample matrix containing the CB count
# create the same table for each class: 
# - all clones
# - co-infection
# - single intection + doublet
# - single infection only

cat("create a clone / sample matrix containing the CB count\n")

# count the clones
clones <- c()
for (i in 1:n_samples) {
    filename <- file.path(sample_dir[i], IN_CLONE_INFO)
    df <- read.table(filename, sep = "\t")
    clones <- c(clones, rownames(df))
}
clones <- unique(clones)

# populate the matrix
M <- matrix(0, nrow = length(clones), ncol = n_samples)
M_single <- matrix(0, nrow = length(clones), ncol = n_samples)
M_single_and_doublet <- matrix(0, nrow = length(clones), ncol = n_samples)
M_coinf <- matrix(0, nrow = length(clones), ncol = n_samples)
rownames(M) <- rownames(M_single) <- rownames(M_single_and_doublet) <- rownames(M_coinf) <- clones
colnames(M) <- colnames(M_single) <- colnames(M_single_and_doublet) <- colnames(M_coinf) <- results_dir_list[,1]
for (i in 1:n_samples) {
    filename <- file.path(sample_dir[i], IN_CLONE_INFO)
    df <- read.table(filename, sep = "\t")
    M[match(rownames(df),rownames(M)),i] <- df$CB.count 
    M_single[match(rownames(df),rownames(M)),i] <- df$CB.count.single
    M_single_and_doublet[match(rownames(df),rownames(M)),i] <- df$CB.count.single + df$CB.count.doublet
    M_coinf[match(rownames(df),rownames(M)),i] <- df$CB.count.coinfected
}

write.table(M, file = OUT_TABLE_CLONE_SAMPLE, sep = "\t", quote = FALSE)
write.table(M_single, file = OUT_TABLE_CLONE_SAMPLE_SINGLE, sep = "\t", quote = FALSE)
write.table(M_single_and_doublet, file = OUT_TABLE_CLONE_SAMPLE_SINGLE_AND_DOUBLET, sep = "\t", quote = FALSE)
write.table(M_coinf, file = OUT_TABLE_CLONE_SAMPLE_COINFECTED, sep = "\t", quote = FALSE)

# generate the barplot with the number of clones
n_clones <- unlist(lapply(1:n_samples, function(x) sum(M[,x] > 0)))
df_clone_bar <- data.frame(sample = sample_name, clone_count = n_clones)
df_clone_bar$sample <- factor(df_clone_bar$sample, levels = sample_name)

write.table(df_clone_bar, file = OUT_TABLE_CLONE_NUM, sep = "\t", row.names = FALSE, quote = FALSE)

g <- ggplot(data = df_clone_bar, aes(x = sample, y = clone_count)) + theme_classic() + geom_bar(stat = "identity") 
g <- g + xlab("") + ylab("Number of clones") + ggtitle("Clones detection")
pdf(OUT_FIGURE_CLONE_NUM, width = 3.5, height = 3)
g
dev.off()

# extract the surviving clones

cat("extract the surviving clones\n")

l <- l_single <- l_single_and_doublet <- l_coinf <- list()
l_freq <- list()
for (i in 1:n_samples) {
    v <- rownames(M)[M[,i] > 0]
    l[[colnames(M)[i]]] <- v
    v <- rownames(M_single)[M_single[,i] > 0]
    l_single[[colnames(M_single)[i]]] <- v
    v <- rownames(M_single_and_doublet)[M_single_and_doublet[,i] > 0]
    l_single_and_doublet[[colnames(M_single_and_doublet)[i]]] <- v
    v <- rownames(M_coinf)[M_coinf[,i] > 0]
    l_coinf[[colnames(M_coinf)[i]]] <- v
    w <- M[M[,i] > 0,i]
    l_freq[[colnames(M)[i]]] <- w
}

venn.diagram(l[idx_surv], filename = OUT_FIGURE_CLONE_SURV_VENN, imagetype = "png")
venn.diagram(l[idx_before_surv], filename = OUT_FIGURE_CLONE_VENN, imagetype = "png")
venn.diagram(l_single[idx_surv], filename = OUT_FIGURE_CLONE_SURV_SINGLE_VENN, imagetype = "png")
venn.diagram(l_single[idx_before_surv], filename = OUT_FIGURE_CLONE_SINGLE_VENN, imagetype = "png")
venn.diagram(l_single_and_doublet[idx_surv], filename = OUT_FIGURE_CLONE_SURV_SINGLE_AND_DOUBLET_VENN, imagetype = "png")
venn.diagram(l_single_and_doublet[idx_before_surv], filename = OUT_FIGURE_CLONE_SINGLE_AND_DOUBLET_VENN, imagetype = "png")
venn.diagram(l_coinf[idx_surv], filename = OUT_FIGURE_CLONE_SURV_COINF_VENN, imagetype = "png")
venn.diagram(l_coinf[idx_before_surv], filename = OUT_FIGURE_CLONE_COINF_VENN, imagetype = "png")

# dotplot of surviving clones frequency
pdf(OUT_FIGURE_CLONE_SURV_CORR, width = 6.5, height = 2.5)
par(mfrow=c(1,length(sample_surv_name)))
for (i in 1:(length(idx_surv)-1)) {
    v1 <- l_freq[[idx_surv[i]]]
    v1_names <- l[[idx_surv[i]]]
    for (j in (i+1):length(idx_surv)) {
        v2 <- l_freq[[idx_surv[j]]]
        v2_names <- l[[idx_surv[j]]]
        cl <- unique(c(v1_names, v2_names))
        x <- y <- rep(0, length(cl))
        x[match(v1_names,cl)] <- v1
        y[match(v2_names,cl)] <- v2
        plot(x, y, xlab = sample_name[idx_surv[i]], ylab = sample_name[idx_surv[j]], cex = 0.5)
    }
}
title("Clone frequency", line = -2, outer = TRUE)
dev.off()

# compute the number of cells for each surviving clones across samples

surv_clones_union <- unique(unlist(l[idx_surv]))
surv_clones_intersection <- surv_clones_union
for (i in idx_surv) {
    surv_clones_intersection <- surv_clones_intersection[surv_clones_intersection %in% l[[i]]]
}

write.table(surv_clones_union, file = OUT_TABLE_SURV_CLONES_UNION, row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(surv_clones_intersection, file = OUT_TABLE_SURV_CLONES_INTERSECTION, row.names = FALSE, col.names = FALSE, quote = FALSE)

# annotate the cells as surviving or not
df_cell_bar <- as.data.frame(matrix(NA, nrow = 0, ncol = 3))
for (i in 1:n_samples) {
    filename <- file.path(sample_dir[i], IN_CELL_INFO)
    df <- read.table(filename, sep = "\t")
    c_surv_union <- c_surv_intersection <- rep("non-surv", nrow(df))
    idx_union <- which(df$expr.GBC.list %in% surv_clones_union)
    c_surv_union[idx_union] <- rep("surv", length(idx_union))
    idx_intersection <- which(df$expr.GBC.list %in% surv_clones_intersection)
    c_surv_intersection[idx_intersection] <- rep("surv", rep = length(idx_intersection))
    df <- cbind(df, surv.union = c_surv_union, surv.intersection = c_surv_intersection)
    write.table(df, file = file.path(sample_dir[i], OUT_CELL_INFO), sep = "\t", quote = FALSE)   
}

# generate the barplot with the classification of clones based on survival after T13
n_cells <- colSums(M)
n_cells
n_surv_clones_union <- colSums(M[clones %in% surv_clones_union,])
n_surv_clones_union
n_surv_clones_intersection <- colSums(M[clones %in% surv_clones_intersection,])
n_surv_clones_intersection
v <- c(n_surv_clones_intersection, n_surv_clones_union - n_surv_clones_intersection, n_cells - n_surv_clones_union)
v
v <- drop(t(as.matrix(v, nrow = length(n_cells), ncol = 3)))
w <- rep(sample_name, 3)
z <- unlist(lapply(c("each T >= 13","any T >= 13","non_surv"), function(x) rep(x, length(sample_name))))
df_cell_surv_bar <- data.frame(sample = w, class = z, cell_count = v)
df_cell_surv_bar$sample <- factor(df_cell_surv_bar$sample, levels = sample_name)
df_cell_surv_bar$class <- factor(df_cell_surv_bar$class, levels = c("each T >= 13","any T >= 13","non_surv"))

write.table(df_cell_surv_bar, file = OUT_TABLE_CELL_SURV_CLASS, sep = "\t", row.names = FALSE, quote = FALSE)

g <- ggplot(data = df_cell_surv_bar, aes(x = sample, y = cell_count, fill = class)) + theme_classic() + geom_bar(stat = "identity") 
g <- g + xlab("") + ylab("Number of cells") + ggtitle("Cell survival", subtitle = "[singlets and doublets with 2 expressed GBC]")
pdf(OUT_FIGURE_CELL_SURV_CLASS, width = 4.5, height = 3)
g
dev.off()

# generate the histogram with the clones grouped by cell numerosity

pdf(OUT_FIGURE_HIST_CLONE_BY_CELL_NUM, width = 3.5*length(idx_before_surv), height = 9)
layout(mat = matrix(1:(3*length(idx_before_surv)), nrow = 3, ncol = length(idx_before_surv)))
for (i in idx_before_surv) {
    c_all <- M[M[,i] > 0,i]
    c_surv_un <- M[M[,i] > 0 & rownames(M) %in% surv_clones_union,i]
    c_surv_int <- M[M[,i] > 0 & rownames(M) %in% surv_clones_intersection,i]
    m <- max(c_all)
    h1 <- hist(c_all, breaks = 0:m, plot = FALSE)
    h2 <- hist(c_surv_un, breaks = 0:m, plot = FALSE)
    h3 <- hist(c_surv_int, breaks = 0:m, plot = FALSE)
    barplot(height = h1$counts, names.arg = 1:m, main = paste(sample_name[i], "\n[all]"), xlab = "number of cells", ylab = "number of clones")
    barplot(height = h2$counts, names.arg = 1:m, main = "[any T >= 13]", xlab = "number of cells", ylab = "number of clones")
    barplot(height = h3$counts, names.arg = 1:m, main = "[each T >= 13]", xlab = "number of cells", ylab = "number of clones")
}
dev.off()

sessionInfo()
q()



