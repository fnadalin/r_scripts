
# R_LIBS_USER="/home/francesca/.local/R-3.6.3/"
# export R_LIBS_USER

# Take a tsv file as input, the column(s) to be added, and an optional prefix for the cell name, and add the info to the object

PREFIX <- ""


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--object", type = "character",
        help="[REQUIRED] input Seurat object"),
    make_option("--table", type = "character",
        help="[REQUIRED] input table in tsv format, with rownames (cells) and colnames (info to be added)"),
    make_option("--columns", type = "character",
        help="[REQUIRED] comma-separated list of field names to be added to the object"),
    make_option("--prefix", type = "character",
        help="[OPTIONAL] prefix to be added to the cell IDs (dash-separated)")
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

if (is.null(opt$table)) {
    write("Option --table is required\nTry --help for help", stderr()) 
    q()
} else {
    TABLE <- opt$table
}

if (is.null(opt$columns)) {
    write("Option --columns is required\nTry --help for help", stderr()) 
    q()
} else {
    COLUMNS <- unlist(strsplit(opt$columns, split = ","))
}

if (!is.null(opt$prefix))
    PREFIX <- opt$prefix


################################## EXECUTION ###################################


library("Seurat")

object <- readRDS(OBJECT)
info <- read.table(TABLE, sep = "\t", header = TRUE)
if (PREFIX != "") 
    rownames(info) <- paste(PREFIX, rownames(info), sep = "-")
if (length(COLUMNS) > 1) { 
    df <- info[,COLUMNS]
} else {
    df <- data.frame(info[,COLUMNS])
    colnames(df) <- COLUMNS
    rownames(df) <- rownames(info)
}
object <- AddMetaData(object, meta = df, col.name = COLUMNS)
saveRDS(object, file = OBJECT)



q()


