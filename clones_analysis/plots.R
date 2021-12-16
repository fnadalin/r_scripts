SUBTITLE <- bquote("[ GBC read fraction" >= " 0.3 ]")
for (s in samples) {
df_cb_gbc_info <- read.table(paste0(s,"/stat_GBC_read_count.tsv"), sep = "\t", header = TRUE)
df_cb_gbc_info <- df_cb_gbc_info[df_cb_gbc_info$rank1.readfrac >= 0.3,]
PLOT_TITLE <- s
OUT_PLOT <- file.path(s, "GBC_read_count_abs_pval_NEW.pdf")
pdf(OUT_PLOT, width = 3.5, height = 4)
plot(x = df_cb_gbc_info$p.val, y = df_cb_gbc_info$readcount, log = "xy", xlab = "p-value", ylab = "GBC read count", main = PLOT_TITLE, sub = SUBTITLE, col = alpha("blue",0.1), pch = 16, cex = 0.5)
abline(v = 1e-10, col = "magenta")
dev.off()
}

MIN_READ_COUNT <- 2
SUBTITLE <- bquote("[ GBC read count" >= ~ .(MIN_READ_COUNT) ~ "]")
for (s in samples) {
df_cb_gbc_info <- read.table(paste0(s,"/stat_GBC_read_count.tsv"), sep = "\t", header = TRUE)
df_cb_gbc_info <- df_cb_gbc_info[df_cb_gbc_info$readcount >= MIN_READ_COUNT,]
PLOT_TITLE <- s
OUT_PLOT <- file.path(s, "GBC_read_count_pval_NEW.pdf")
pdf(OUT_PLOT, width = 3.5, height = 4)
plot(x = df_cb_gbc_info$p.val, y = df_cb_gbc_info$rank1.readfrac, log = "x", xlab = "p-value", ylab = "(GBC read count) / (rank1 GBC read count)", main = PLOT_TITLE, sub = SUBTITLE, col = alpha("blue",0.1), pch = 16, cex = 0.5)
abline(v = 1e-20, col = "magenta")
abline(h = 0.3, col = "green")
dev.off()
}

