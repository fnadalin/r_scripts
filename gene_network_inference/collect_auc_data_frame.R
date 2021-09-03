
# create the regulon / cell matrix from separate AUCell runs (on the same SCENIC input!!!)
# entries of the matrix are the AUC values
# each input should contain the same regulons, in the same order

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript collect_auc_data_frame.R <scenic_folder_list> <out_rds>\n\n")
	cat("<scenic_folder_list>   file with the names of SCENIC folders (containing int/), one per line\n")
	cat("<out_rds>              output object containing the regulon / cell data frame\n\n")
	q()
}

SCENIC_IN_DIR <- args[1]
OUT_RDS <- args[2]

# load the libraries
library("SCENIC")
library("AUCell")

scenic_in_dir <- drop(as.matrix(read.table(SCENIC_IN_DIR)))

df <- data.frame()
for (dir in scenic_in_dir) {
	aucell_regulonAUC <- readRDS(file.path(dir, "int", "3.4_regulonAUC.Rds"))
	auc <- getAUC(aucell_regulonAUC)
	if (ncol(df) == 0) {
		df <- as.data.frame(auc)
	} else {
		df <- cbind(df, as.data.frame(auc))
	}
}

saveRDS(df, file = OUT_RDS)

q()

