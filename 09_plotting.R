#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("2 arguments are needed:\n",
      "  1) list of smooth.co.txt for tested population (one path per line)\n",
      "  2) file with IDs to filter out (blacklist)\n", sep = "")
  stop("Requires command line arguments.")
}

suppressMessages(library(dplyr))

# -------------------------
# Configuration parameters
# -------------------------

CHR_ENDS <- c(30427671, 19698289, 23459830, 18585056, 26975502)

BIN_SIZE <- 200000
SMOOTH_K <- 5

`%notin%` <- function(x, y) !(x %in% y)

process_filenames <- function(filename) {
  first_part_removed  <- sub("^[^._]*[._]", "", filename)
  second_part_removed <- sub("\\..*$", "", first_part_removed)
  second_part_removed
}

# -------------------------
# Input reading
# -------------------------

list1 <- read.table(args[1], stringsAsFactors = FALSE)  # tested: paths (or filenames)
filt  <- na.omit(read.table(args[2], fill = TRUE))

# -------------------------
# Parse refined.corrected.breaks files into CO events
# -------------------------

parse_smooth_file <- function(file_path, lib_id, n_chr = 5) {
  dat <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)

  all <- NULL
  for (i in 1:n_chr) {
    chr <- dat[dat[, 1] == i, , drop = FALSE]
    if (nrow(chr) > 1) {
      start <- chr[, 2]
      stop  <- chr[, 3]

      stop  <- stop[-1]
      start <- start[-length(start)]

      width <- stop - start
      cos   <- start + round(width / 2)

      lib  <- rep(lib_id, length(cos))
      chrs <- rep(i, length(cos))

      bind <- cbind(lib, chrs, start, stop, cos, width)
      all <- rbind(all, bind)
    }
  }
  all
}

# Tested population
taf4b.cos <- NULL
for (k in 1:nrow(list1)) {
  file_path <- as.character(list1[k, 1])
  if (!file.exists(file_path)) stop(paste("File does not exist:", file_path))
  lib_id <- as.character(list1[k, 1])  # keep the identifier as written in the list
  taf4b.cos <- rbind(taf4b.cos, parse_smooth_file(file_path, lib_id, n_chr = 5))
}

# Add dummy librarynum, convert, clean lib IDs
taf4b.cos <- cbind(librarynum = rep(1, nrow(taf4b.cos)), taf4b.cos)
taf4b.cos <- data.frame(taf4b.cos, stringsAsFactors = FALSE)
taf4b.cos$lib <- as.character(taf4b.cos$lib)
# Enforce expected column names
colnames(taf4b.cos) <- c("librarynum","lib","chrs","start","stop","cos","width")

# Numeric columns
taf4b.cos$chrs  <- as.numeric(taf4b.cos$chrs)
taf4b.cos$start <- as.numeric(taf4b.cos$start)
taf4b.cos$stop  <- as.numeric(taf4b.cos$stop)
taf4b.cos$cos   <- as.numeric(taf4b.cos$cos)
taf4b.cos$width <- as.numeric(taf4b.cos$width)

# -------------------------
# Filter individuals and emit Histogram.txt + COs.txt
# -------------------------

count_per_lib <- function(cos_df) {
  uni <- unique(as.character(cos_df$lib))
  all.ind <- vapply(uni, function(id) {
    nrow(cos_df[cos_df$lib == id, , drop = FALSE])
  }, numeric(1))
  data.frame(COs = as.numeric(all.ind), ID = uni, stringsAsFactors = FALSE)
}

lol <- count_per_lib(taf4b.cos)
lol <- lol[which(lol$COs < (mean(lol$COs) + 2 * sd(lol$COs))), ]
lol <- lol[which(lol$COs > (mean(lol$COs) - 2 * sd(lol$COs))), ]
#lol <- lol[lol$ID %notin% filt$V1, ]

write.table(lol, "Histogram.txt", row.names = FALSE, sep = "\t", quote = FALSE)

taf4b.cos <- taf4b.cos[taf4b.cos$lib %in% lol$ID, ]
write.table(taf4b.cos, "COs.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Denominator for smoothed counts
n_taf4b <- length(unique(taf4b.cos$lib))
if (n_taf4b == 0) stop("No tested libraries found after filtering (n_taf4b=0).")
cat("[info] Tested libraries:", n_taf4b, "\n")

# -------------------------
# Across-chromosomes profile + emission
# -------------------------

ma <- rep(1 / SMOOTH_K, SMOOTH_K)

fill_edge_na <- function(x) {
  x <- as.numeric(x)
  if (!any(is.na(x))) return(x)
  ok <- which(!is.na(x))
  if (length(ok) == 0) return(x)
  x[1:(min(ok) - 1)] <- x[min(ok)]
  x[(max(ok) + 1):length(x)] <- x[max(ok)]
  x
}

tha.cum <- c(0, cumsum(CHR_ENDS))

cos.matrix <- NULL
for (i in 1:5) {
  taf4b.chr <- taf4b.cos[taf4b.cos$chrs == i, , drop = FALSE]

  windows <- seq(1, CHR_ENDS[i], by = BIN_SIZE)
  if (tail(windows, 1) != CHR_ENDS[i]) windows <- c(windows, CHR_ENDS[i])
  starts <- windows[-length(windows)]
  ends   <- windows[-1]

  cum.windows <- starts + tha.cum[i]

  taf_pos_col <- "cos"
  taf4b_pos <- as.numeric(taf4b.chr[[taf_pos_col]])

  win.taf4b <- numeric(length(starts))
  for (j in 1:length(starts)) {
    win.taf4b[j] <- sum(is.finite(taf4b_pos) & taf4b_pos >= starts[j] & taf4b_pos <= ends[j])
  }
  win.taf4b <- win.taf4b / n_taf4b
  filt.taf4b <- fill_edge_na(stats::filter(win.taf4b, ma))

  chr.dat <- cbind(rep(i, length(cum.windows)), cum.windows, filt.taf4b)
  cos.matrix <- rbind(cos.matrix, chr.dat)
}

cos.matrix2 <- cos.matrix
colnames(cos.matrix2) <- c("chr", "cum.windows", "filt.taf4b")

write.table(cos.matrix2,
            "Across_chromosomes.txt",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

# -------------------------
# Plot across-chromosomes profile
# -------------------------

plot_across_chromosomes <- function(cos_mat, chr_ends,
                                    out_png = "Across_chromosomes.png",
                                    width_px = 1600, height_px = 700, res = 150) {
  cos_mat <- as.data.frame(cos_mat)

  x  <- suppressWarnings(as.numeric(cos_mat[, 2]))
  y  <- suppressWarnings(as.numeric(cos_mat[, 3]))

  y_fin <- y[is.finite(y)]
  if (length(y_fin) == 0) stop("No finite y-values to plot.")

  y_max <- max(y_fin)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1

  chr_cum_ends <- cumsum(chr_ends)
  chr_cum_starts <- c(0, chr_cum_ends[-length(chr_cum_ends)])
  chr_mid <- (chr_cum_starts + chr_cum_ends) / 2
  chr_lab <- paste0("Chr", seq_along(chr_ends))

  png(out_png, width = width_px, height = height_px, res = res)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mar = c(5, 5, 3, 2))

  plot(x, y, type = "l",
       xlab = "", ylab = "CO rate per bin (per individual, smoothed)",
       ylim = c(0, y_max * 1.05), xaxt = "n")

  abline(v = chr_cum_ends, lty = 3)
  axis(1, at = chr_mid, labels = chr_lab, tick = FALSE)

  invisible(TRUE)
}

plot_across_chromosomes(cos.matrix2, CHR_ENDS, out_png = "Across_chromosomes.png")
