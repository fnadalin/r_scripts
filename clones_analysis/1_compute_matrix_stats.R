
# Input:
#  - sparse CRISPR matrix
#  - min read count for GBC detection 
# Output:
#  - table with CB, GBC, read.count information (for each (CB,GBC) pair)
#  - plot with the cumulative read frequency across GBC
#  - plot with the GBC read count and GBC read count / CB read count

MIN_READ_COUNT <- 5
PLOT_TITLE <- ""
MIN_UMI_COUNT <- 0

FEATURE_ASSAY <- "CRISPR Guide Capture"
FEATURE_NAME <- "GBC"

################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--matrix_dir", type = "character",
        help="[REQUIRED] directory containing matrix.mtx, features.tsv, barcodes.tsv, where features are both genes and GBC (may be gzipped)"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] output folder containing tables and plots"),
    make_option("--min_UMI_count", type = "integer", default = MIN_UMI_COUNT,
        help="[REQUIRED] minimum UMI count to consider a CB as detected [default=%default]"),
    make_option("--min_read_count", type = "integer", default = MIN_READ_COUNT,
        help="[OPTIONAL] minimum number of reads to consider a GBC as detected [default=%default]"),
    make_option("--plot_title", type = "character",
        help="[OPTIONAL] title to be written on the plots"),
    make_option("--feature_assay", type = "character",
        help="[OPTIONAL] name of the feature assay (\"CRISPR Guide Capture\", \"MULTISEQ Capture\"...)")	
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$matrix_dir)) {
    write("Option --matrix_dir is required\nTry --help for help", stderr()) 
    q()
} else {
    MATRIX_DIR <- opt$matrix_dir
}

if (is.null(opt$out_dir)) {
    write("Option --out_dir is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_DIR <- opt$out_dir
}

if (!is.null(opt$min_UMI_count))
    MIN_UMI_COUNT <- opt$min_UMI_count

if (!is.null(opt$min_read_count))
    MIN_READ_COUNT <- opt$min_read_count

if (!is.null(opt$plot_title))
    PLOT_TITLE <- opt$plot_title

if (!is.null(opt$feature_assay)) 
    FEATURE_ASSAY <- opt$feature_assay

if (FEATURE_ASSAY == "MULTISEQ Capture")
    FEATURE_NAME <- "MBC"

################################### IMPORT #####################################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "../functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

# print parameters to file
sink(file = file.path(OUT_DIR, paste0("parameters_", SCRIPT_NAME)))
for (i in 1:length(opt))
	cat(paste0(names(opt)[i], "=", opt[[i]], "\n"))
sink()


################################## EXECUTION ###################################

# output file names

if (!dir.exists(OUT_DIR))
    dir.create(OUT_DIR, recursive = TRUE)

OUT_TABLE_CB_STATS        <- file.path(OUT_DIR, "stat_CB_UMI_read_count.tsv")
OUT_TABLE_GBC_STATS       <- file.path(OUT_DIR, paste0("stat_",FEATURE_NAME,"_read_count.tsv"))
OUT_TABLE_GBC             <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_read_count.tsv"))
OUT_FIGURE_CB_UMI         <- file.path(OUT_DIR, "CB_UMI_count_density.pdf")
OUT_FIGURE_GBC            <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_read_count_cum.pdf"))
OUT_FIGURE_GBC_READ_COUNT_PVAL <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_read_count_pval.pdf"))
OUT_FIGURE_GBC_READ_COUNT_DENSITY <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_read_count_density.pdf"))
OUT_FIGURE_GBC_PVAL_DENSITY <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_pval_density.pdf"))
OUT_FIGURE_GBC_PVAL_HIST <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_pval_hist.pdf"))
OUT_FIGURE_GBC_READFRAC_DENSITY <- file.path(OUT_DIR, paste0(FEATURE_NAME,"_readfrac_density.pdf"))
OUT_FIGURE_CB_GBC         <- file.path(OUT_DIR, paste0("CB_",FEATURE_NAME,"_read_count.pdf"))
OUT_FIGURE_CB_GBC_UMI     <- file.path(OUT_DIR, paste0("CB_",FEATURE_NAME,"_read_UMI_count.pdf"))
OUT_FIGURE_CB_GBC_READ_COUNT_DENSITY <- file.path(OUT_DIR, paste0("CB_",FEATURE_NAME,"_read_count_density.pdf"))
OUT_FIGURE_CB_RANK1_GBC_READ_COUNT_DENSITY <- file.path(OUT_DIR, paste0("CB_rank1_",FEATURE_NAME,"_read_count_density.pdf"))

# parse expression matrix

cat("parse input files\n")

M_genes <- Matrix()
M_gbc <- Matrix()
# ParseFeatureBarcodeMatrixCombined(MATRIX_DIR, M_genes, M_gbc)
ParseFeatureBarcodeMatrixExtract(MATRIX_DIR, M_genes, "Gene Expression")
ParseFeatureBarcodeMatrixExtract(MATRIX_DIR, M_gbc, FEATURE_ASSAY)
cb_umicount <- colSums(M_genes)
M_genes <- M_genes[,cb_umicount > MIN_UMI_COUNT]
M_gbc <- M_gbc[,cb_umicount > MIN_UMI_COUNT]
cb_umicount <- cb_umicount[cb_umicount > MIN_UMI_COUNT]

cb_exp <- colnames(M_genes)
genes <- rownames(M_genes)
# cb_genecount <- colSums(M_genes > 0) # NEW: count the number of detected genes
cb_gbc <- colnames(M_gbc)
gbc <- rownames(M_gbc)

if (length(cb_exp) != length(cb_gbc)) {
    write("Gene expression and GBC expression matrices must contain the same number of CB", stderr()) 
    q()
}

compare <- unlist(lapply(1:length(cb_exp), function(x) cb_gbc[x] != cb_exp[x]))
if (sum(compare) > 1) {
    write("Gene expression and GBC expression matrices must contain the same CB, in the same order", stderr()) 
    q()
}

cb <- cb_exp

# count the number of reads for each (CB,GBC) pair, only for those satisfying read.count >= MIN_READ_COUNT

cat("count the number of reads for each (CB,GBC) pair, only for those satisfying read.count >= MIN_READ_COUNT\n")

s <- summary(M_gbc)
s <- s[s[,3] >= MIN_READ_COUNT,]
readcount <- s[,3]
CB <- cb_gbc[s[,2]]
GBC <- gbc[s[,1]]
num_reads <- sum(readcount)

# compute the top ranked guide for each cell

cat("compute the top ranked guide for each cell\n")

idx <- order(readcount, decreasing = TRUE)
readcount <- readcount[idx]
CB <- CB[idx]
GBC <- GBC[idx]

rank1_idx <- match(cb, CB)
cb_rank1_idx <- match(CB[rank1_idx], cb)
cb_idx <- match(CB, cb)

df_cb_gbc_info <- data.frame(CB = CB, GBC = GBC, readcount = readcount)
rownames(df_cb_gbc_info) <- unlist(lapply(1:length(CB), function(x) paste(CB[x], GBC[x], sep = ".")))
df_cb_gbc_info <- df_cb_gbc_info[order(df_cb_gbc_info$CB),]

# compute the total UMI count and the total read count (from guides with read.count >= MIN_READ_COUNT) across the cells

cat("compute the total UMI count and the total read count (from guides with read.count >= MIN_READ_COUNT) across the cells\n")

cb_readcount <- unlist(lapply(cb, function(x) sum(readcount[CB == x])))
cb_rank1_gbc <- GBC[rank1_idx]
cb_rank1_gbc <- cb_rank1_gbc[cb_rank1_idx]
cb_rank1_gbc_readcount <- readcount[rank1_idx]
cb_rank1_gbc_readcount <- cb_rank1_gbc_readcount[cb_rank1_idx]
df_cb_info <- data.frame(UMIcount = cb_umicount, readcount = cb_readcount, rank1.GBC = cb_rank1_gbc, rank1.GBC.readcount = cb_rank1_gbc_readcount)
rownames(df_cb_info) <- cb

write.table(df_cb_info, file = OUT_TABLE_CB_STATS, quote = FALSE, sep = "\t")

g <- ggplot(data = df_cb_info, aes(x = UMIcount)) + theme_classic() + geom_density() + xlab("UMI count") + ylab("density") + ggtitle(PLOT_TITLE)
pdf(OUT_FIGURE_CB_UMI, width = 3.5, height = 3.5)
g
dev.off()

pdf(OUT_FIGURE_CB_GBC, width = 5.5, height = 6)
SUBTITLE <- bquote(paste0("[ ",FEATURE_NAME," read count") >= ~ .(MIN_READ_COUNT) ~ "]")
plot(x = readcount, y = readcount / cb_readcount[cb_idx], log = "x",
     xlab = "GBC read count", ylab = paste0("(",FEATURE_NAME," read count) / (CB read count)"), main = PLOT_TITLE, sub = SUBTITLE,
     col = "black", pch = 16, cex = 0.5)
points(x = readcount[rank1_idx], y = readcount[rank1_idx] / cb_readcount[cb_rank1_idx],
       col = "magenta", pch = 16, cex = 0.5)
legend("bottomright", legend = c("rank = 1", "rank >= 1"), col = c("magenta", "black"), pch = 16, pt.cex = 0.5, cex = 0.8)
dev.off()

pdf(OUT_FIGURE_CB_GBC_UMI, width = 5.5, height = 6)
SUBTITLE <- bquote(paste0("[ ",FEATURE_NAME," read count") >= ~ .(MIN_READ_COUNT) ~ "]")
plot(x = readcount, y = cb_umicount[cb_idx], log = "xy",
     xlab = paste0(FEATURE_NAME," read count"), ylab = "CB UMI count", main = PLOT_TITLE, sub = SUBTITLE,
     col = "black", pch = 16, cex = 0.5)
points(x = readcount[rank1_idx], y = cb_umicount[cb_rank1_idx],
       col = "magenta", pch = 16, cex = 0.5)
legend("bottomright", legend = c("rank = 1", "rank >= 1"), col = c("magenta", "black"), pch = 16, pt.cex = 0.5, cex = 0.8)
dev.off()

# compute the p-value of a guide to be expressed

cat("compute the p-value of a guide to be expressed\n")

# extract the info of the fraction of expressed reads in a guide (wrt to the cell)
idx <- match(df_cb_gbc_info$CB, rownames(df_cb_info))
gbc_readfrac <- df_cb_gbc_info$readcount / df_cb_info$readcount[idx]
rank1_gbc_readfrac <- df_cb_gbc_info$readcount / df_cb_info$rank1.GBC.readcount[idx]

# Poisson test
pval <- rep(0, nrow(df_cb_gbc_info))
ncells <- nrow(df_cb_info)
gbc_unique <- unique(df_cb_gbc_info$GBC)
gbc_mean_read_count <- unlist(lapply(gbc_unique, function(x) {
                              idx <- which(df_cb_gbc_info$GBC == x)
                              sum(df_cb_gbc_info$readcount[idx]/ncells) 
                              }))
names(gbc_mean_read_count) <- gbc_unique
pval <- unlist(lapply(1:nrow(df_cb_gbc_info), function(x) {
               idx <- which(names(gbc_mean_read_count) == df_cb_gbc_info$GBC[x])
               pois <- poisson.test(df_cb_gbc_info$readcount[x], r = gbc_mean_read_count[idx], alternative = "g")
               pois$p.value
               }))
names(pval) <- df_cb_gbc_info$GBC

df_cb_gbc_info <- cbind(df_cb_gbc_info, data.frame(readfrac = gbc_readfrac, rank1.readfrac = rank1_gbc_readfrac, p.val = pval))

write.table(df_cb_gbc_info, file = OUT_TABLE_GBC_STATS, quote = FALSE, sep = "\t")

pdf(OUT_FIGURE_GBC_READ_COUNT_PVAL, width = 3.5, height = 3.5)
SUBTITLE <- bquote(paste0("[ ",FEATURE_NAME," read count") >= ~ .(MIN_READ_COUNT) ~ "]")
plot(x = df_cb_gbc_info$p.val, y = df_cb_gbc_info$rank1.readfrac, log = "x",
     xlab = "p-value", ylab = paste0("(",FEATURE_NAME," read count) / (rank1 ",FEATURE_NAME," read count)", main = PLOT_TITLE, sub = SUBTITLE,
     col = "black", pch = 16, cex = 0.5)
dev.off()

g <- ggplot(data = df_cb_gbc_info, aes(x = readcount)) + theme_classic() + geom_density() + xlab(paste0(FEATURE_NAME," read count")) + ylab("density") + ggtitle(PLOT_TITLE)
pdf(OUT_FIGURE_GBC_READ_COUNT_DENSITY, width = 3.5, height = 3.5)
g
dev.off()

df_cb_gbc_info$p.val[df_cb_gbc_info$p.val == 0] <- rep(1e-308, sum(df_cb_gbc_info$p.val == 0))
g <- ggplot(data = df_cb_gbc_info, aes(x = -log10(p.val))) + theme_classic() + geom_density() + xlab("-log10(p-value)") + ylab("density") + ggtitle(PLOT_TITLE)
pdf(OUT_FIGURE_GBC_PVAL_DENSITY, width = 3.5, height = 3.5)
g
dev.off()

# generate the histogram as an alternative to the density
pdf(OUT_FIGURE_GBC_PVAL_HIST, width = 3.5, height = 3.5)
m <- max(-log10(df_cb_gbc_info$p.val))
h <- hist(-log10(df_cb_gbc_info$p.val), breaks = 51, plot = FALSE)
barplot(height = h$counts, names.arg = h$breaks[2:length(h$breaks)], las = 2, main = PLOT_TITLE, xlab = "-log10(p-value)", ylab = paste0("number of ",FEATURE_NAME))
dev.off()

g <- ggplot(data = df_cb_gbc_info, aes(x = rank1.readfrac)) + theme_classic() + geom_density() + xlab(paste0("(",FEATURE_NAME," read count) / (rank1 ",FEATURE_NAME," read count)")) + ylab("density") + ggtitle(PLOT_TITLE)
pdf(OUT_FIGURE_GBC_READFRAC_DENSITY, width = 3.5, height = 3.5)
g
dev.off()

# compute the total read count across the guides

cat("compute the total read count across the guides\n")

gbc_readcount <- unlist(lapply(gbc, function(x) sum(readcount[GBC == x])))
names(gbc_readcount) <- gbc
gbc_readcount <- gbc_readcount[gbc_readcount > 0]
gbc_readcount <- gbc_readcount[order(gbc_readcount, decreasing = TRUE)]

write.table(gbc_readcount, file = OUT_TABLE_GBC, sep = "\t", quote = FALSE, col.names = FALSE)

# compute the cumulative distribution of the total read count across the guides

cat("compute the cumulative distribution of the total read count across the guides\n")

num_detected_gbc <- sum(gbc_readcount > 0)
cum_gbc_readcount <- rep(0, num_detected_gbc+1)
for (i in 1:num_detected_gbc) {
    cum_gbc_readcount[i+1] <- cum_gbc_readcount[i] + gbc_readcount[i]
}
cum_gbc_readcount <- cum_gbc_readcount[2:length(cum_gbc_readcount)]

pdf(OUT_FIGURE_GBC, width = 4, height = 4.5)
SUBTITLE <- bquote(paste0("[ ",FEATURE_NAME," read count") >= ~ .(MIN_READ_COUNT) ~ "]")
plot(x = 1:length(cum_gbc_readcount), y = cum_gbc_readcount / num_reads, 
     lty = 1, col = "blue", 
     main = PLOT_TITLE, sub = SUBTITLE, xlab = paste0("Number of ",FEATURE_NAME), ylab = "Fraction of reads")
dev.off() 

g <- ggplot(data = df_cb_info, aes(x = readcount)) + theme_classic() + geom_density() + xlab(paste0("CB ",FEATURE_NAME," read count")) + ylab("density") + ggtitle(PLOT_TITLE)
pdf(OUT_FIGURE_CB_GBC_READ_COUNT_DENSITY, width = 3.5, height = 3.5)
g
dev.off()

g <- ggplot(data = df_cb_info, aes(x = rank1.GBC.readcount)) + theme_classic() + geom_density() + xlab(paste0("rank1 ",FEATURE_NAME," read count")) + ylab("density") + ggtitle(PLOT_TITLE)
pdf(OUT_FIGURE_CB_RANK1_GBC_READ_COUNT_DENSITY, width = 3.5, height = 3.5)
g
dev.off()


sessionInfo()
q()



