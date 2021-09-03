# Extract the names of the genes for each regulon

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript regulon_info.R <scenic_folder> <out_file>]\n\n")
	cat("<scenic_folder>     SCENIC folder (containing int/ and output/)\n")
	cat("<out_file>          output .tsv file where the list of genes are printed\n\n")
	q()
}

# parse the input parameters
SCENIC_FOLDER <- args[1]
OUT_FILE <- args[2]

reg <- readRDS(file.path(SCENIC_FOLDER, "int", "3.1_regulons_forAUCell.Rds"))
x <- unlist(lapply(names(reg), function(x) paste(unlist(reg[x]), collapse=",")))

df <- data.frame(regulon = names(reg), genes = x)

write.table(df, file = OUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)

q()

