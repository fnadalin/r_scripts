
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript generate_random_cell_subset.R <in_cells> <seed> <num_cells> <out_cells>\n\n")
	cat("<in_cells>      list of cell IDs, one per line\n")
	cat("<seed>          for the rng\n")
	cat("<num_cells>     cells to generate\n")
	cat("<out_cells>     list of cell IDs\n\n")
	q()
}

# parse the input parameters
IN_FILE <- args[1]
SEED <- as.numeric(args[2])
NUM_CELLS <- as.numeric(args[3])
OUT_FILE <- args[4]

cells <- drop(as.matrix(read.table(IN_FILE)))

set.seed(SEED)
s <- sample(x = 1:length(cells), size = NUM_CELLS)
s_cells <- cells[sort(s)]

write.table(s_cells, file = OUT_FILE, quote = FALSE, col.names = FALSE, row.names = FALSE)

sessionInfo()
q()

