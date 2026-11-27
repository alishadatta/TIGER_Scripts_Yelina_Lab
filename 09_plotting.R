#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

cat("Starting TIGER Recombination Analysis from HMM output (rough_COs_ERR*.txt)...\n")

# 1. Target rough_COs_ERR*.txt files (Column 6 holds the Viterbi HMM state)
base_dir <- "/scratch/ltwb0962/Tiger_Pipeline/tiger_input_output/tiger500"

file_list <- list.files(
  path = base_dir,
  pattern = "^rough_COs_ERR[0-9]+\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(file_list) == 0) {
  stop("No rough_COs_ERR*.txt files found! Check base_dir.")
}

cat(paste("Found", length(file_list), "F2 HMM output files to process.\n"))

# 2. Extract HMM blocks, remove short single-marker glitches (< 150 kb), and call crossovers
all_co_events <- data.frame()
valid_samples <- 0
min_hmm_block_bp <- 150000

for (f in file_list) {
  tryCatch({
    if (file.exists(f) && file.info(f)$size > 0) {
      df <- read.table(f, header = FALSE, sep = "", stringsAsFactors = FALSE)

      if (nrow(df) >= 2 && ncol(df) >= 6) {
        df_hmm <- data.frame(
          Chr = suppressWarnings(as.numeric(df[, 2])),
          Pos = suppressWarnings(as.numeric(df[, 3])),
          Genotype = as.character(df[, 6]),
          stringsAsFactors = FALSE
        ) %>%
          filter(!is.na(Chr), !is.na(Pos), Chr %in% 1:5, Genotype %in% c("CC", "CL", "LL")) %>%
          arrange(Chr, Pos)

        if (nrow(df_hmm) > 0) {
          valid_samples <- valid_samples + 1

          # Collapse consecutive markers with the same HMM genotype into blocks
          df_blocks <- df_hmm %>%
            group_by(Chr) %>%
            mutate(block_id = cumsum(Genotype != lag(Genotype, default = "START"))) %>%
            group_by(Chr, block_id, Genotype) %>%
            summarise(Start = min(Pos), End = max(Pos), .groups = "drop") %>%
            # Remove short HMM glitches (< 150 kb)
            filter((End - Start) >= min_hmm_block_bp)

          # Re-merge adjacent blocks that now share the same genotype
          if (nrow(df_blocks) >= 2) {
            df_blocks <- df_blocks %>%
              group_by(Chr) %>%
              mutate(clean_id = cumsum(Genotype != lag(Genotype, default = "START"))) %>%
              group_by(Chr, clean_id, Genotype) %>%
              summarise(Start = min(Start), End = max(End), .groups = "drop") %>%
              arrange(Chr, Start)

            co_idx <- which(
              df_blocks$Chr[-nrow(df_blocks)] == df_blocks$Chr[-1] &
              df_blocks$Genotype[-nrow(df_blocks)] != df_blocks$Genotype[-1]
            )

            if (length(co_idx) > 0) {
              co_events <- data.frame(
                Chr = df_blocks$Chr[co_idx],
                CO_Position = (df_blocks$End[co_idx] + df_blocks$Start[co_idx + 1]) / 2
              )
              all_co_events <- rbind(all_co_events, co_events)
            }
          }
        }
      }
    }
  }, error = function(e) {})
}

if (nrow(all_co_events) == 0) {
  stop("No crossover events extracted.")
}

cat(paste("Valid F2 samples:", valid_samples,
          "| Total COs:", nrow(all_co_events),
          "| Mean COs/F2:", round(nrow(all_co_events) / valid_samples, 2),
          "(Target: ~7.80)\n"))

# 3. Bin into fixed 200-kb windows (0.2 Mb)
window_size_bp <- 200000
window_size_mb <- 0.2

all_co_events$Window_Mb <- floor(all_co_events$CO_Position / window_size_bp) * window_size_mb + (window_size_mb / 2)

binned_data <- all_co_events %>%
  group_by(Chr, Window_Mb) %>%
  summarise(Count = n(), .groups = "drop")

chr_lengths <- c(30.4, 19.7, 23.5, 18.5, 27.0)
full_grid <- data.frame()
for (ch in 1:5) {
  max_w <- ceiling(chr_lengths[ch] / window_size_mb) * window_size_mb
  windows <- seq(window_size_mb / 2, max_w, by = window_size_mb)
  full_grid <- rbind(full_grid, data.frame(Chr = ch, Window_Mb = windows))
}

binned_data <- full_grid %>%
  left_join(binned_data, by = c("Chr", "Window_Mb")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
  mutate(Count = case_when(
    Chr == 4 & Window_Mb >= 1.8 & Window_Mb <= 4.8 ~ 0, # Mask Ler heterochromatic inversion
    TRUE ~ Count
  ))

# Divide by 2 * valid_samples gametes
binned_data <- binned_data %>%
  mutate(CO_Rate = (Count / (2 * valid_samples)) * 100 / window_size_mb)

binned_data$Chr <- factor(binned_data$Chr, levels = 1:5, labels = paste("Chr", 1:5))

# 4. Generate the Final Publication Plot
p <- ggplot(binned_data, aes(x = Window_Mb, y = CO_Rate)) +
  geom_line(color = "#2b7bba", linewidth = 0.6) +
  facet_grid(Chr ~ .) +
  scale_x_continuous(limits = c(0, 31), breaks = seq(0, 30, 10)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(
    x = "Position (Mb)",
    y = "CO rate (cM/Mb)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#e6e6e6", color = "black", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.spacing = unit(0.5, "lines")
  )

output_file <- "population_crossover_landscape_HMM.pdf"
ggsave(output_file, plot = p, width = 9, height = 10, dpi = 300)
cat(paste("Success! HMM-based plot saved as:", output_file, "\n"))
