#!/usr/bin/env Rscript

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: plotting.R <file_list.txt>")
}

files <- readLines(args[1])
files <- files[file.exists(files)]

cat("Break files found:", length(files), "\n")

CHR_ENDS <- c(
  30427671,
  19698289,
  23459830,
  18585056,
  26975502
)

CENTROMERES <- data.frame(
  chr = 1:5,
  start = c(14000000, 3000000, 13000000, 3500000, 11000000),
  end   = c(16000000, 5000000, 15000000, 5500000, 13000000)
)

BUFFER <- 10000

parse_breaks <- function(file) {

  dat <- tryCatch(
    read.table(
      file,
      header = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) return(NULL)
  )

  if (is.null(dat) || ncol(dat) < 4)
    return(NULL)

  dat <- dat[,1:4]

  colnames(dat) <- c(
    "chr",
    "start",
    "end",
    "state"
  )

  dat$chr <- suppressWarnings(as.numeric(dat$chr))
  dat$start <- suppressWarnings(as.numeric(dat$start))
  dat$end <- suppressWarnings(as.numeric(dat$end))

  dat <- dat[complete.cases(dat), ]

  dat <- dat[
    dat$state == "CL" &
    dat$end > dat$start,
  ]

  if (nrow(dat) == 0)
    return(NULL)

  dat$pos <- (dat$start + dat$end) / 2

  return(dat)
}

all_cos <- do.call(
  rbind,
  lapply(files, parse_breaks)
)

if (is.null(all_cos) || nrow(all_cos) == 0) {
  stop("No CL intervals found.")
}

cat(
  "CL intervals loaded:",
  nrow(all_cos),
  "\n"
)

for (i in 1:nrow(CENTROMERES)) {

  chr <- CENTROMERES$chr[i]
  s <- CENTROMERES$start[i] - BUFFER
  e <- CENTROMERES$end[i] + BUFFER

  all_cos <- all_cos[!(
    all_cos$chr == chr &
    all_cos$pos >= s &
    all_cos$pos <= e
  ), ]
}

cat(
  "After centromere filtering:",
  nrow(all_cos),
  "\n"
)

offset <- c(
  0,
  cumsum(CHR_ENDS)
)[1:5]

all_cos$genome_pos <-
  all_cos$pos +
  offset[all_cos$chr]

dens <- density(
  all_cos$genome_pos,
  bw = 100000
)

png(
  "CL_density_genomewide.png",
  width = 1800,
  height = 700,
  res = 150
)

plot(
  dens,
  main = "Genome-wide Distribution of Recombinant (CL) Intervals",
  xlab = "Genome position",
  ylab = "Density",
  lwd = 2
)

abline(
  v = cumsum(CHR_ENDS),
  lty = 3,
  col = "grey60"
)

axis(
  1,
  at = cumsum(CHR_ENDS) - CHR_ENDS/2,
  labels = paste0("Chr", 1:5),
  tick = FALSE,
  line = 1
)

dev.off()

cat(
  "Plot written: CL_density_genomewide.png\n"
)
