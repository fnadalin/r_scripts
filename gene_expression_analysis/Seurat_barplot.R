# Created on 2020/09/15

########################### DEFAULT PARAMETER VALUES ###########################

DO_SUBSET <- FALSE

################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--object", type = "character",
        help="[REQUIRED] .Robj file containing the Seurat object"),
    make_option("--out_prefix", type = "character",
        help="[REQUIRED] output file containing the plot"),
    make_option("--cond_mode", type = "character",
        help="[REQUIRED] mode ID for barplot conditions"),
    make_option("--split_mode", type = "character",
        help="[REQUIRED] mode ID for additional conditions for splitting the bars"),
    make_option("--subset", type = "character",
        help="[OPTIONAL] identity mode to subset by"),
    make_option("--subset_id", type = "character",
        help="[OPTIONAL] ID for the subset")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$object)) {
    write("Option --object is required\nTry --help for help", stderr()) 
    q()
} else {
    OBJECT <- opt$object
}

if (is.null(opt$out_prefix)) {
    write("Option --out_prefix is required\nTry --help for help", stderr()) 
    q()
} else {
    OUT_PREFIX <- opt$out_prefix
}

if (is.null(opt$cond_mode)) {
    write("Option --cond_mode is required\nTry --help for help", stderr()) 
    q()
} else {
    COND_MODE <- opt$cond_mode
}

if (is.null(opt$split_mode)) {
    write("Option --split_mode is required\nTry --help for help", stderr()) 
    q()
} else {
    SPLIT_MODE <- opt$split_mode
}

DO_SUBSET <- !(is.null(opt$subset) | is.null(opt$subset_id))
if (DO_SUBSET) {
    SUBSET <- opt$subset
    SUBSET_ID <- opt$subset_id
}


################################## EXECUTION ###################################

library("ggplot2")

object <- readRDS(OBJECT)

data <- data.frame(object@meta.data[[COND_MODE]], object@meta.data[[SPLIT_MODE]])
colnames(data) <- c(COND_MODE, SPLIT_MODE)
if (DO_SUBSET) {
    idx <- which(object@meta.data[[SUBSET]] == SUBSET_ID)
    data <- data[idx,]
}

cond <- sort(unique(data[[COND_MODE]]))
split <- sort(unique(data[[SPLIT_MODE]]))

df <- matrix(NA, nrow = length(cond)*length(split), ncol = 3)
df <- as.data.frame(df)
colnames(df) <- c("cond", "split", "num")
df$num <- 0 

n <- 1
for (s in cond) {
    for (c in split) {
        count <- sum(data[[COND_MODE]] == s & data[[SPLIT_MODE]] == c)
        df[n,] <- c(s, c, count)
        n <- n + 1
    }
}
    
df$num <- as.numeric(df$num)
write.table(df, file = paste0(OUT_PREFIX, ".tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

g <- ggplot(df, aes(x=cond, y=num, fill=split)) + theme_minimal() + xlab("") + theme_classic() 
if (DO_SUBSET)
    g <- g + ggtitle(SUBSET_ID)
g <- g + geom_bar(stat="identity", position=position_dodge()) + ylab("number of cells")
pdf(paste0(OUT_PREFIX, ".pdf"), width = 4, height = 3)
print(g)
dev.off()



q()


