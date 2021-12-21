
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
        help="[REQUIRED] tsv file containing the files with the CB classification for each sample, one per line, where each column represents:\n		1 - the sample name\n		2 - the directory containing the results of the CB classification"),
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

library("ggplot2")

if (!dir.exists(OUT_DIR))
    dir.create(OUT_DIR, recursive = TRUE)

IN_CLONE_INFO              <- "stats_clones.tsv"
IN_CELL_INFO               <- "CB_classification.tsv"

OUT_TABLE_CLONE_SAMPLE     <- file.path(OUT_DIR, "clone_sample_CB_count.tsv")
OUT_TABLE_CLONE_SAMPLE_SINGLE <- file.path(OUT_DIR, "clone_sample_CB_count_single.tsv")
OUT_TABLE_CLONE_SAMPLE_SINGLE_AND_DOUBLET <- file.path(OUT_DIR, "clone_sample_CB_count_single_and_doublet.tsv")
OUT_TABLE_CLONE_SAMPLE_COINFECTED <- file.path(OUT_DIR, "clone_sample_CB_count_coinfected.tsv")
OUT_TABLE_CELL_CLASS       <- file.path(OUT_DIR, "barplot_CB_classification.tsv")
OUT_TABLE_CLONE_NUM        <- file.path(OUT_DIR, "barplot_GBC_count.tsv")
OUT_FIGURE_CELL_CLASS      <- file.path(OUT_DIR, "barplot_CB_classification.pdf")
OUT_FIGURE_CELL_CLASS_FRAC <- file.path(OUT_DIR, "barplot_CB_classification_frac.pdf")
OUT_FIGURE_CLONE_NUM       <- file.path(OUT_DIR, "barplot_GBC_count.pdf")

# parse the input

results_dir_list <- as.matrix(read.table(RESULTS_DIR_LIST, sep = "\t"))
sample_name <- results_dir_list[,1]
sample_dir <- results_dir_list[,2]

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
g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(OUT_FIGURE_CELL_CLASS, width = 4.5, height = 3)
g
dev.off()

g <- ggplot(data = df_cell_bar, aes(x = sample, y = count, fill = class)) + theme_classic() + geom_bar(stat = "identity", position = "fill") 
g <- g + xlab("") + ylab("Number of CB") + ggtitle("CB classification")
g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
g <- g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(OUT_FIGURE_CLONE_NUM, width = 3.5, height = 3)
g
dev.off()


sessionInfo()
q()



