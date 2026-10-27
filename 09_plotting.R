#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

cat("Starting TIGER Recombination Analysis (Sliding-Window Block Parser)...\n")

# 1. Dynamically target the sliding_window.breaks.txt files
base_dir <- "/scratch/ltwb0962/Tiger_Pipeline/input/tiger500"

err_folders <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)
err_folders <- err_folders[grep("ERR[0-9]+$", err_folders)]

file_list <- c()
for(folder in err_folders) {
  target_files <- list.files(path = folder,
                             pattern = "INDEX_.*_sliding_window\\.breaks\\.txt$",
                             full.names = TRUE)
  file_list <- c(file_list, target_files)
}

if (length(file_list) == 0) {
  stop("No sliding_window.breaks.txt files found! Double check path structures.")
}

cat(paste("Found", length(file_list), "valid genotype block files to process.\n"))

# 2. Parse transitions/crossovers between genotype states
all_co_events <- data.frame()
total_samples <- length(file_list)

for (f in file_list) {
  tryCatch({
    if (file.exists(f) && file.info(f)$size > 0) {
      df <- read.table(f, header = FALSE, sep = "", stringsAsFactors = FALSE)

      if (nrow(df) >= 2 && ncol(df) >= 5) {
        # Columns: 1=Sample, 2=Chr, 3=Start, 4=End, 5=Genotype
        df_clean <- df[, 2:5]
        colnames(df_clean) <- c("Chr", "Start", "End", "Genotype")

        df_clean$Chr <- as.numeric(df_clean$Chr)
        df_clean$Start <- as.numeric(df_clean$Start)
        df_clean$End <- as.numeric(df_clean$End)

        # Sort sequentially along chromosomes
        df_clean <- df_clean %>% arrange(Chr, Start)

        co_events <- data.frame()
        for(i in seq_len(nrow(df_clean) - 1)) {
          # Crossover occurs if it's the same chromosome but the genotype state switches
          if(df_clean$Chr[i] == df_clean$Chr[i+1] && df_clean$Genotype[i] != df_clean$Genotype[i+1]) {
            # Crossover position is the junction point between blocks
            co_pos <- df_clean$End[i]

            co_events <- rbind(co_events, data.frame(
              Chr = df_clean$Chr[i],
              CO_Position = co_pos
            ))
          }
        }

        if(nrow(co_events) > 0) {
          all_co_events <- rbind(all_co_events, co_events)
        }
      }
    }
  }, error = function(e) {
    # Ignore anomalous lines
  })
}

if (nrow(all_co_events) == 0) {
  stop("Aggregation stopped: No state transitions found in sliding window data.")
}

cat(paste("Successfully isolated", nrow(all_co_events), "true population crossover junctions.\n"))

# 3. Bin into fixed 200-kb windows (0.2 Mb) matching the figure legend
window_size_bp <- 200000
window_size_mb <- 0.2

all_co_events$Window_Mb <- floor(all_co_events$CO_Position / window_size_bp) * window_size_mb + (window_size_mb / 2)

binned_data <- all_co_events %>%
  group_by(Chr, Window_Mb) %>%
  summarise(Count = n(), .groups = 'drop')

# Complete the grid for zero-crossover windows across Arabidopsis chromosomes
chr_lengths <- c(30.4, 19.7, 23.5, 18.5, 27.0)
full_grid <- data.frame()
for(ch in 1:5) {
  max_w <- ceiling(chr_lengths[ch] / window_size_mb) * window_size_mb
  windows <- seq(window_size_mb/2, max_w, by = window_size_mb)
  full_grid <- rbind(full_grid, data.frame(Chr = ch, Window_Mb = windows))
}

binned_data <- full_grid %>%
  left_join(binned_data, by = c("Chr", "Window_Mb")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count))

# Calculate accurate CO rate metric (cM/Mb)
binned_data <- binned_data %>%
  mutate(CO_Rate = (Count / total_samples) * 100 / window_size_mb)

binned_data$Chr <- factor(binned_data$Chr, levels = 1:5, labels = paste("Chr", 1:5))

# 4. Generate the Final Publication Plot
p <- ggplot(binned_data, aes(x = Window_Mb, y = CO_Rate)) +
  geom_line(color = "#2b7bba", linewidth = 0.6) +
  facet_grid(Chr ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(0, 31), breaks = seq(0, 30, 10)) +
  labs(
    x = "Position (Mb)",
    y = "CO rate (cM/Mb)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "#e6e6e6", color = "grey80"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93"),
    panel.spacing = unit(0.5, "lines")
  )

# 5. Save Output
output_file <- "population_crossover_landscape.pdf"
ggsave(output_file, plot = p, width = 9, height = 10, dpi = 300)
cat(paste("Success! Clean publication-ready plot saved as:", output_file, "\n"))
