
# Input:
#  - table with (GBC,CB) info on read count
#  - table with CB UMI count and read count info
#  - min/max UMI count to retain a CB
#  to classify a guide as expressed:
#  - min number of reads per guide
#  - min GBC read count / CB read count
#  - whether to use or not a Poisson distribution  
# Output:
#  - augmented table for (CB,GBC) pairs with info on whether the GBC is expressed in CB or not
#  - table with CB, expressed GBC list, quality pass (based on UMI count), type (unclassified, uninf, single inf, coinf, doublet)
#  - plot with the read count / UMI count across CB, coloured by class

FEATURE_ASSAY <- "CRISPR Guide Capture"
FEATURE_NAME <- "GBC"

MIN_UMI_COUNT <- 0
MAX_UMI_COUNT <- 100000
MIN_READ_COUNT_EXPR <- 20
MIN_READ_FRAC_EXPR <- 0.2
PLOT_TITLE <- ""
MIN_RANK1_READ_FRAC_EXPR <- 0.3
PVAL_CUTOFF <- 1e-5
MIN_CLONE_FREQ_COINFECTED <- 5 # minimum number of cells where the clone is detected to safely classify it as such in a co-infection event
DOUBLET_PROB_CUTOFF <- 0.05


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
    make_option("--doublet_rate", type = "double",
        help="[REQUIRED] expected doublet rate, given the number of loaded cells"),
    make_option("--min_umi_count", type = "integer", default = MIN_UMI_COUNT,
        help="[OPTIONAL] minimum UMI count for quality pass [default=%default]"),
    make_option("--max_umi_count", type = "integer", default = MAX_UMI_COUNT,
        help="[OPTIONAL] maximum UMI count for non-doublet cells [default=%default]"),
    make_option("--min_read_count_expr", type = "integer", default = MIN_READ_COUNT_EXPR,
        help="[OPTIONAL] minimum read count to classify a GBC as expressed [default=%default]"),
    make_option("--min_read_frac_expr", type = "double", default = MIN_READ_FRAC_EXPR,
        help="[OPTIONAL] minimum GBC read count / CB read count to classify a GBC as expressed [default=%default]"),
    make_option("--min_rank1_read_frac_expr", type = "double", default = MIN_RANK1_READ_FRAC_EXPR,
        help="[OPTIONAL] minimum GBC read count / rank 1 GBC read count to classify a GBC as expressed [default=%default]"),
    make_option("--pval_cutoff", type = "double", default = PVAL_CUTOFF,
        help="[OPTIONAL] maximum p-value of the Poisson test to accept a GBC as expressed [default=%default]"),
    make_option("--doublet_prob_cutoff", type = "double", default = DOUBLET_PROB_CUTOFF,
        help="[OPTIONAL] maximum probability for the doublet test to label a CB as doublet [default=%default]"),
    make_option("--plot_title", type = "character",
        help="[OPTIONAL] title to be written on the plots"),
    make_option("--feature_assay", type = "character",
        help="[OPTIONAL] name of the feature assay (\"CRISPR Guide Capture\", \"MULTISEQ Capture\"...)")

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

if (is.null(opt$doublet_rate)) {
    write("Option --doublet_rate is required\nTry --help for help", stderr()) 
    q()
} else {
    DOUBLET_RATE <- opt$doublet_rate
}

if (!is.null(opt$min_umi_count))
    MIN_UMI_COUNT <- opt$min_umi_count

if (!is.null(opt$max_umi_count))
    MAX_UMI_COUNT <- opt$max_umi_count

if (!is.null(opt$min_read_count_expr))
    MIN_READ_COUNT_EXPR <- opt$min_read_count_expr

if (!is.null(opt$min_read_frac_expr))
    MIN_READ_FRAC_EXPR <- opt$min_read_frac_expr

if (!is.null(opt$min_rank1_read_frac_expr))
    MIN_RANK1_READ_FRAC_EXPR <- opt$min_rank1_read_frac_expr

if (!is.null(opt$pval_cutoff))
    PVAL_CUTOFF <- opt$pval_cutoff

if (!is.null(opt$doublet_prob_cutoff))
    DOUBLET_PROB_CUTOFF <- opt$doublet_prob_cutoff

if (!is.null(opt$plot_title))
    PLOT_TITLE <- opt$plot_title

if (!is.null(opt$feature_assay))
    FEATURE_ASSAY <- opt$feature_assay

if (FEATURE_ASSAY == "MULTISEQ Capture")
    FEATURE_NAME <- "MBC"

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

library("quantreg")
library("ggplot2")
library("scales")

# output file names

if (!dir.exists(OUT_DIR))
    dir.create(OUT_DIR, recursive = TRUE)

IN_CB_GBC_INFO_TABLE <- file.path(IN_DIR, paste0("stat_",FEATURE_NAME,"_read_count.tsv"))
IN_CB_INFO_TABLE     <- file.path(IN_DIR, "stat_CB_UMI_read_count.tsv")

OUT_TABLE_GBC_STATS  <- file.path(OUT_DIR, paste0("stat_",FEATURE_NAME,"_read_count_with_info.tsv"))
OUT_TABLE_STATS      <- file.path(OUT_DIR, "CB_classification.tsv")
OUT_FIGURE_DOUBLET_SCORE <- file.path(OUT_DIR, "CB_multi_guide_doublet_score_density.pdf")
OUT_FIGURE_CB_GBC    <- file.path(OUT_DIR, "CB_classification.pdf")

# parse input files

df_cb_gbc_info <- read.table(IN_CB_GBC_INFO_TABLE, sep = "\t")
df_cb_info <- read.table(IN_CB_INFO_TABLE, sep = "\t")

is_expressed <- (df_cb_gbc_info$p.val < PVAL_CUTOFF & df_cb_gbc_info$readfrac >= MIN_READ_FRAC_EXPR & df_cb_gbc_info$rank1.readfrac >= MIN_RANK1_READ_FRAC_EXPR & df_cb_gbc_info$readcount >= MIN_READ_COUNT_EXPR)
df_cb_gbc_info <- cbind(df_cb_gbc_info, is.expressed = is_expressed)

write.table(df_cb_gbc_info, file = OUT_TABLE_GBC_STATS, sep = "\t", quote = FALSE)

# detect the expressed guides for each CB

cat("detect the expressed guides for each CB\n")

# extract the list of expressed guides for each CB
cb_with_expr_guides <- unique(df_cb_gbc_info$CB[is_expressed])
num_expressed_gbc <- unlist(lapply(cb_with_expr_guides, function(x) sum(is_expressed & df_cb_gbc_info$CB == x)))
expressed_gbc <- unlist(lapply(cb_with_expr_guides, function(x) {
                        paste(sort(df_cb_gbc_info$GBC[is_expressed & df_cb_gbc_info$CB == x]), collapse = ",")
                        }))

cb_gbc_list <- rep("", nrow(df_cb_info))
cb_gbc_list[match(cb_with_expr_guides, rownames(df_cb_info))] <- expressed_gbc

cb_gbc_num <- rep(0, nrow(df_cb_info))
cb_gbc_num[match(cb_with_expr_guides, rownames(df_cb_info))] <- num_expressed_gbc

# classify the CB 

cat("classify the CB\n")

# CB quality check 
quality_pass <- (df_cb_info$UMIcount >= MIN_UMI_COUNT)

idx_single_guide <- which(cb_gbc_num == 1)
idx_doublet_high_umi <- which(df_cb_info$UMIcount > MAX_UMI_COUNT)

# look at the read count / UMI count ratio to classify the cells as doublets or co-infection events
N <- nrow(df_cb_info)
D <- ceiling(DOUBLET_RATE*N/100) # total num of doublets
h <- length(idx_doublet_high_umi) # cells > MAX_UMI_COUNTS are classified as doublets
d <- D - h # remaining doublets to call

cat(paste("all cells:",    nrow(df_cb_info), "\n"))
cat(paste("quality pass:", sum(quality_pass), "\n"))
cat(paste("single guide:", sum(cb_gbc_num == 1), "\n"))
cat(paste("no exp guide:", sum(cb_gbc_num == 0), "\n"))
cat(paste("multi guide:",  sum(cb_gbc_num > 1), "\n"))
cat(paste("doublets (> MAX_UMI_COUNT):", h, "\n"))
cat(paste("doublets (<= MAX_UMI_COUNT):", d, "\n"))

n <- sum(cb_gbc_num > 1 & df_cb_info$UMIcount <= MAX_UMI_COUNT) # number of multi-guide cells
r <- d/n # doublet_rate
cat(paste("doublet rate across multi-guide cells:", r, "\n"))
cat(paste("expect", d, "doublets overall\n"))

# quantile regression on the multi-guide cells
# top: co-infected, bottom: doublet
# fit the cells with m expressed guides independently 
idx_coinfected <- c()
idx_doublet_low_ratio <- c()
prob_doublet <- rep(0, length(cb_gbc_list))
for (m in 2:max(cb_gbc_num)) {
    cat(paste("analyze the set of", m, "guides\n"))
    idx_multi_guide <- which(cb_gbc_num == m & df_cb_info$UMIcount <= MAX_UMI_COUNT)
    nn <- length(idx_multi_guide)
    cat(paste(nn, "cells found\n"))
    cc <- floor(nn*(1-r)) # number of co-infected cells
    cat(paste("expect", cc, "co-infected cells among the cells with", m, "expressed guides\n"))
    dd <- floor(nn*r) # number of doublets
    cat(paste("expect", dd, "doublets among the cells with", m, "expressed guides\n"))
    v <- c()
    if (length(idx_multi_guide) > 1) {
        # quantile_fit <- rq(readcount ~ UMIcount, data = df_cb_info[idx_multi_guide,], tau = r)
        # idx_coinfected <- c(idx_coinfected, idx_multi_guide[quantile_fit$residuals > 0]) # OLD
        # idx_doublet_low_ratio <- c(idx_doublet_low_ratio, idx_multi_guide[quantile_fit$residuals < 0]) # OLD
        # starting from the top residual, iteratively include new clones (and the relative cells) corresponding to the current set of guides 
        x <- cb_gbc_list[idx_multi_guide]
        # sort the GBC sets to select the top read count / UMI count cells as co-infected clones
        x <- unique(x[order(df_cb_info$readcount[idx_multi_guide] / df_cb_info$UMIcount[idx_multi_guide], decreasing = TRUE)]) 
        c <- 0 # cumulative number of cells for each GBC set
        classified_cells <- length(c(idx_single_guide, idx_coinfected))
        for (i in 1:length(x)) {
            clone_freq <- sum(cb_gbc_list == x[i])
            prob <- 0
#            if (clone_freq < MIN_CLONE_FREQ_COINFECTED) {
                # generate all possible pairs of clones from the guides (1- and 2-guides clones are enough)
                # compute the frequency of the two clones in the single infection and in the already computed co-infection events
                cl <- unlist(strsplit(x[i], split = ","))
                for (j in 1:length(cl)) {
                    c1 <- cl[j]
                    c2 <- paste(cl[!(1:length(cl) == j)], collapse = ",")
                    f1 <- sum(cb_gbc_list[c(idx_single_guide, idx_coinfected)] == c1)
                    f2 <- sum(cb_gbc_list[c(idx_single_guide, idx_coinfected)] == c2)
                    pr_part <- (f1+f2)/classified_cells
                    if (length(cl) == 2)
                        pr_part <- pr_part / 2
                    prob <- prob + pr_part
                }
                if (length(cl) > 2) {
                    for (j in 1:(length(cl)-1)) {
                        for (k in (i+1):length(cl)) {
                            c1 <- paste(cl[c(j,k)], collapse = ",")
                            c2 <- paste(cl[!(1:length(cl) %in% c(j,k))], collapse = ",")
                            f1 <- sum(cb_gbc_list[c(idx_single_guide, idx_coinfected)] == c1)
                            f2 <- sum(cb_gbc_list[c(idx_single_guide, idx_coinfected)] == c2)
                            pr_part <- (f1+f2)/classified_cells
                            if (length(cl) == 4)
                                pr_part <- pr_part / 2
                            prob <- prob + pr_part
                        }
                    }
                }
#            }
            prob_doublet[which(cb_gbc_list == x[i])] <- min(prob,1)
            if ((clone_freq >= MIN_CLONE_FREQ_COINFECTED) || (prob < DOUBLET_PROB_CUTOFF)) {
                c <- c + clone_freq
                # label each cell with this GBC set as co-infected
                v <- c(v, which(cb_gbc_list == x[i]))
                if (c >= cc) 
                    break
            }
        }
    } 
    cat(paste(length(v), "cells found as co-infected based on top read count / UMI count sets of guides\n"))
    cat(paste(sum(!(idx_multi_guide %in% v)), "cells found as doublets\n"))
    idx_coinfected <- c(idx_coinfected, v)
    idx_doublet_low_ratio <- c(idx_doublet_low_ratio, idx_multi_guide[!(idx_multi_guide %in% v)])
}

# plot the density of the doublet score
doublet_score = prob_doublet[cb_gbc_num > 1 & df_cb_info$UMIcount <= MAX_UMI_COUNT]
df <- data.frame(doublet_score = doublet_score)
q <- quantile(doublet_score, probs = 0.99)
cat(paste("99th quantile of the doublet score:", q, "\n"))
q <- scientific(q, digits = 2)
SUBTITLE <-  paste("[ 99th percentile =", q, "]")
g <- ggplot(data = df, aes(x = doublet_score)) + theme_classic() + geom_density() + ggtitle(PLOT_TITLE, sub = SUBTITLE)
pdf(OUT_FIGURE_DOUBLET_SCORE, width = 3.5, height = 3.5)
g
dev.off()

idx_doublet <- c(idx_doublet_low_ratio, idx_doublet_high_umi)

class <- rep("uninfected", nrow(df_cb_info))
class[idx_single_guide] <- rep("single_guide", length(idx_single_guide))
class[idx_coinfected] <- rep("coinfected", length(idx_coinfected))
class[idx_doublet] <- rep("doublet", length(idx_doublet))

df <- data.frame(UMIcount = df_cb_info$UMIcount, readcount = df_cb_info$readcount, quality.pass = quality_pass, 
                 expr.GBC.num = cb_gbc_num, expr.GBC.list = cb_gbc_list, prob.doublet = prob_doublet, class = class)
rownames(df) <- rownames(df_cb_info)

write.table(df, file = OUT_TABLE_STATS, sep = "\t", quote = FALSE)

pdf(OUT_FIGURE_CB_GBC, width = 5, height = 5.5)
plot(x = df_cb_info$UMIcount, y = df_cb_info$readcount, col = "gray", pch = 16, cex = 0.3, main = PLOT_TITLE, xlab = "UMI count", ylab = "read count")
points(x = df_cb_info$UMIcount[idx_single_guide], y = df_cb_info$readcount[idx_single_guide], col = "blue", pch = 16, cex = 0.3)
points(x = df_cb_info$UMIcount[idx_coinfected], y = df_cb_info$readcount[idx_coinfected], col = "red", pch = 16, cex = 0.3)
points(x = df_cb_info$UMIcount[idx_doublet], y = df_cb_info$readcount[idx_doublet], col = "green", pch = 16, cex = 0.3)
legend("topright", legend = c("uninfected", "single-infection", "co-infection", "doublet"), pch = c(16,16,16,16),
       col = c("gray", "blue", "red", "green"), lty = 0, bty = "n")
dev.off()



sessionInfo()
q()



