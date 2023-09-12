# rename duplicated features

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <in> <out>\n")
    cat("(10x output format)\n")
    q()
}

IN <- args[1]
OUT <- args[2]

if (length(grep("gz",IN)) > 0) {
    features <- read.table(gzfile(IN), sep = "\t")
} else {
    features <- read.table(IN, sep = "\t")
}

with_header <- FALSE
if (length(grep("ENS",features[1,1])) == 0) {
    with_header <- TRUE
    colnames(features) <- features[1,]
    features <- features[2:nrow(features),]
}

uniq_features <- unique(features[,2])
uniq_features_freq <- unlist(lapply(uniq_features, function(x) sum(features[,2] == x)))
idx <- which(uniq_features_freq > 1)
uniq_features <- uniq_features[idx]
uniq_features_freq <- uniq_features_freq[idx]
for (x in 1:length(uniq_features)) {
    idx <- which(features[,2] == uniq_features[x])
    features[idx,2] <- paste(uniq_features[x],1:uniq_features_freq[x],sep=".")
}

write.table(features, file = OUT, sep = "\t", quote = FALSE, row.names = FALSE, col.names = with_header)
if (length(grep("gz",OUT)) > 0) {
    system(paste0("gzip", OUT))
}

q()

