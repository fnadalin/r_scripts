
# Input:
#  - table with CB, expressed GBC list, quality pass (based on UMI count), type (unclassified, uninf, single inf, coinf, doublet)
# Output:
#  - table with the list of clones and the number of cells per clone
#  - plot with the cumulative number of cells per clone

PLOT_TITLE <- ""


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--in_dir", type = "character",
        help="[REQUIRED] input folder containing pre-computed stats"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] output folder containing tables and plots"),
    make_option("--plot_title", type = "character",
        help="[OPTIONAL] title to be written on the plots")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$in_dir)) {
    write("Option --in_dir is required\nTry --help for help", stderr()) 
    q()
} else {
    IN_DIR <- opt$in_dir
}

if (is.null(opt$out_dir)) {
    write("Option --out_dir is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_DIR <- opt$out_dir
}

if (!is.null(opt$plot_title))
    PLOT_TITLE <- opt$plot_title


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

if (!dir.exists(OUT_DIR))
    dir.create(OUT_DIR, recursive = TRUE)

IN_CB_INFO_TABLE       <- file.path(IN_DIR, "CB_classification.tsv")

OUT_TABLE_CLONES_STATS <- file.path(OUT_DIR, "stats_clones.tsv")
OUT_FIGURE_CLONES_CB   <- file.path(OUT_DIR, "clones_CB_count.pdf")

# parse the input file

df_cb_info <- read.table(IN_CB_INFO_TABLE, sep = "\t")
df_cb_info$expr.GBC.list <- as.character(df_cb_info$expr.GBC.list)
df_cb_info$expr.GBC.num <- as.numeric(df_cb_info$expr.GBC.num) 

# extract the clones

cat("extract the clones\n")

# frequencies of clones from single-infections
idx_single <- which(df_cb_info$expr.GBC.num > 0 & df_cb_info$class == "single_guide")
v <- df_cb_info$expr.GBC.list[idx_single]
freq <- as.matrix(table(v))
cl_single <- rownames(freq)
cl_freq_single <- c(freq)
names(cl_freq_single) <- cl_single

# frequencies of clones from co-infections
idx_co <- which(df_cb_info$expr.GBC.num > 0 & df_cb_info$class == "coinfected")
v <- df_cb_info$expr.GBC.list[idx_co]
freq <- as.matrix(table(v))
cl_co <- rownames(freq)
cl_freq_co <- c(freq)
names(cl_freq_co) <- cl_co

# frequencies of clones from doublets (<= 2 GBC, otherwise assignment of CB to clones is not unique)
idx_doublet <- which(df_cb_info$expr.GBC.num > 0 & df_cb_info$expr.GBC.num <= 2 & df_cb_info$class == "doublet")
s <- paste(df_cb_info$expr.GBC.list[idx_doublet], collapse = ",")
v <- unlist(strsplit(s, split = ","))
freq <- as.matrix(table(v))
cl_doublet <- rownames(freq)
cl_freq_doublet <- c(freq)
names(cl_freq_doublet) <- cl_doublet

# collect all clones
v <- c(cl_freq_single, cl_freq_co, cl_freq_doublet)
cl_all <- unique(names(v))
n <- length(cl_all)

# collect the frequencies per clone, divided by CB class
M <- matrix(0, nrow = n, ncol = 4)
df <- as.data.frame(M)
colnames(df) <- c("CB.count", "CB.count.single", "CB.count.coinfected", "CB.count.doublet")
rownames(df) <- cl_all
df$CB.count.single[match(cl_single, cl_all)] <- cl_freq_single
df$CB.count.coinfected[match(cl_co, cl_all)] <- cl_freq_co
df$CB.count.doublet[match(cl_doublet, cl_all)] <- cl_freq_doublet
df$CB.count <- rowSums(df)

# print the CB count per clone

df <- df[order(df$CB.count, decreasing = TRUE),]
write.table(df, file = OUT_TABLE_CLONES_STATS, sep = "\t", quote = FALSE)

# compute the cumulative distribution of the CB count across the clones

cat("compute the cumulative distribution of the CB count across the clones\n")

num_cells <- nrow(df_cb_info)
cum_clone_cb_count <- rep(0, n+1)
for (i in 1:n) {
    cum_clone_cb_count[i+1] <- cum_clone_cb_count[i] + df$CB.count[i]
}
cum_clone_cb_count <- cum_clone_cb_count[2:length(cum_clone_cb_count)]

pdf(OUT_FIGURE_CLONES_CB, width = 4, height = 4.5)
SUBTITLE <- cat(paste("from", num_cells, "cells\n"))
plot(x = 1:length(cum_clone_cb_count), y = cum_clone_cb_count / num_cells, 
     lty = 1, col = "blue", 
     main = PLOT_TITLE, sub = SUBTITLE, xlab = "Number of clones", ylab = "Fraction of CB")
dev.off() 



sessionInfo()
q()



