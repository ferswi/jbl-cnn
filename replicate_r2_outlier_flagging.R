#!/usr/bin/env Rscript
# ============================================================================
# Flag R^2 outlier peaks in the replicate concordance analysis.
#
# For each replicate pair, a peak's distance from the y=x line is d = log10(a)-log10(b).
# Peaks whose robust z-score |(d - median)/MAD| exceeds `nmad` are flagged as outliers
# (they sit off the main cloud in the R^2 scatter). Peaks flagged in MULTIPLE pairs
# are "consistent" outliers -- the recurring off-line pattern -- and are ranked highest.
#
# Reuses the replicate_qc outputs (signal/*.tab + consensus_windows.bed). Base R only.
#
# Usage:
#   Rscript replicate_r2_outlier_flagging.R <signal_dir> <windows_bed> <outdir> \
#           <celltype> <nmad> rep1 rep2 ...
#   e.g.  ... replicate_qc/ex_endo_parietal/signal \
#             replicate_qc/ex_endo_parietal/consensus_windows.bed \
#             replicate_qc/ex_endo_parietal ex_endo_parietal 3.5 1a 1b 2a 2b
# ============================================================================

args     <- commandArgs(trailingOnly = TRUE)
sig_dir  <- args[1]
win_bed  <- args[2]
outdir   <- args[3]
celltype <- args[4]
nmad     <- as.numeric(args[5])
reps     <- args[-(1:5)]

# --- region x replicate matrix from bigWigAverageOverBed tabs (mean0 = col 5) ----
tabs <- lapply(reps, function(r) {
  d <- read.table(file.path(sig_dir, paste0(r, ".tab")),
                  header = FALSE, stringsAsFactors = FALSE)[, c(1, 5)]
  names(d) <- c("region", r); d
})
m <- Reduce(function(a, b) merge(a, b, by = "region"), tabs)
rownames(m) <- m$region
M  <- as.matrix(m[, reps, drop = FALSE])
pc <- 0.01
L  <- log10(M + pc)

# --- genomic coordinates from consensus_windows.bed (col4 = region name) ---------
win <- read.table(win_bed, header = FALSE, stringsAsFactors = FALSE)
colnames(win)[1:4] <- c("chrom", "start", "end", "region")
coord <- win[match(rownames(M), win$region), c("chrom", "start", "end")]

# --- per-pair outlier detection --------------------------------------------------
n          <- length(reps)
flag_count <- setNames(integer(nrow(M)), rownames(M))  # in how many pairs is a peak an outlier
maxdist    <- numeric(nrow(M))                          # largest |residual| across pairs
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  d      <- L[, i] - L[, j]
  s      <- mad(d)
  z      <- if (s > 0) (d - median(d)) / s else rep(0, length(d))
  is_out <- abs(z) > nmad
  flag_count <- flag_count + as.integer(is_out)
  maxdist    <- pmax(maxdist, abs(d))
}

# --- outputs ---------------------------------------------------------------------
npairs <- n * (n - 1) / 2
out <- data.frame(coord, region = rownames(M),
                  n_pairs_outlier = flag_count,
                  max_log10_diff  = round(maxdist, 3),
                  round(M, 4), check.names = FALSE)
out <- out[out$n_pairs_outlier > 0, ]
out <- out[order(-out$n_pairs_outlier, -out$max_log10_diff), ]
write.csv(out, file.path(outdir, paste0(celltype, "_r2_outliers.csv")), row.names = FALSE)

# consistent outliers = flagged in at least half the pairs -> a BED for inspection
consistent <- out[out$n_pairs_outlier >= ceiling(npairs / 2),
                  c("chrom", "start", "end", "region", "n_pairs_outlier")]
write.table(consistent, file.path(outdir, paste0(celltype, "_r2_outliers_consistent.bed")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# scatter grid with outliers highlighted in red
pdf(file.path(outdir, paste0(celltype, "_r2_outliers_scatter.pdf")), width = 8, height = 8)
panel.pts <- function(x, y, ...) {
  d <- x - y; s <- mad(d)
  z <- if (s > 0) (d - median(d)) / s else rep(0, length(d))
  col <- ifelse(abs(z) > nmad, rgb(1, 0, 0, 0.6), rgb(0, 0, 0, 0.2))
  points(x, y, pch = ".", col = col); abline(0, 1, col = "blue", lty = 2)
}
pairs(L, labels = reps, panel = panel.pts, gap = 0.3,
      main = paste0(celltype, " - R2 outliers (red = >", nmad, " MAD off y=x)"))
dev.off()

# --- console summary -------------------------------------------------------------
cat("consensus regions:        ", nrow(M), "\n")
cat("threshold:                ", nmad, "MAD from the median residual\n")
cat("regions flagged in >=1 pair:", nrow(out), "\n")
cat("consistent (>= half of", npairs, "pairs):", nrow(consistent), "\n")
cat("outputs written to", outdir, ":\n")
cat("  ", paste0(celltype, "_r2_outliers.csv"), "  (all flagged, ranked)\n")
cat("  ", paste0(celltype, "_r2_outliers_consistent.bed"), "  (recurring outliers)\n")
cat("  ", paste0(celltype, "_r2_outliers_scatter.pdf"), "  (scatter, outliers in red)\n")
