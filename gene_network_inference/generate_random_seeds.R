
MAX <- 1000000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript generate_random_seeds.R <rng> <num_seeds> <out_file>\n\n")
	cat("<rng>         seed for the random number generator\n")
	cat("<num_seeds>   number of random seeds to generate\n")
	cat("<out_file>    where the random seeds are saved\n\n")
	q()
}

# parse the input parameters
RNG <- as.numeric(args[1])
NUM_SEEDS <- as.numeric(args[2])
OUT_FILE <- args[3]

set.seed(RNG)
s <- sample(x = 1:MAX, size = NUM_SEEDS)

write.table(s, file = OUT_FILE, quote = FALSE, row.names = FALSE, col.names = FALSE)

sessionInfo()
q()

