#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

# -----------------------------
# STEP 1: Summarize recombination/deletion data
# -----------------------------


recomb_files <- list.files(pattern = "sorted\\.txt$")

process_file <- function(file) {
  message("Summarizing: ", file)
  data <- read_tsv(file, show_col_types = FALSE)
  
  summarized_data <- data %>%
    group_by(chrom, start, end) %>%
    summarize(
      depth              = n(),
      depth_positive     = sum(strand == "+"),
      depth_negative     = sum(strand == "-"),
      n_fragments        = n_distinct(read_name),
      max_left_overhang  = max(left_overhang),
      max_right_overhang = max(right_overhang),
      max_aligned_length = max(aligned_length),
      .groups = "drop"
    )
  
  # Examples:
  #   filtered_..._deletion.sorted.txt       -> filtered_..._deletion.count.txt
  #   filtered_..._recombinations.sorted.txt -> filtered_..._recombinations.count.txt
  output_file <- sub("sorted\\.txt$", "count.txt", file)
  write_tsv(summarized_data, file = output_file)
  message("  -> wrote: ", output_file)
}

invisible(lapply(recomb_files, process_file))

# -----------------------------
# STEP 2: Frequency calculation and filtering
# -----------------------------

coverage_files <- list.files(path = ".", pattern = "coverage$", full.names = TRUE)
count_files    <- list.files(path = ".", pattern = "count\\.txt$", full.names = FALSE)

for (coverage_path in coverage_files) {
  cov_base <- basename(coverage_path)
  message("\nProcessing coverage: ", cov_base)
  
  count_file <- NULL
  
  # ---------- CASE 1: STAR-style ----------

  if (grepl("annotated_standardizationAligned\\.out\\.sorted\\.coverage$", cov_base)) {
    # Define a "core" string to match in the count files
    core <- sub("annotated_standardizationAligned\\.out\\.sorted\\.coverage$",
                "annotated_standardizationAligned", cov_base)
    # Look for any count file containing this core and "deletion"
    candidates <- count_files[grepl(core, count_files) & grepl("deletion", count_files)]
    if (length(candidates) > 0) {
      count_file <- candidates[1]
    }
  }
  
  # ---------- CASE 2: ViReMa-style ----------
  
  if (is.null(count_file) && grepl("recombinations\\.coverage$", cov_base)) {
    core <- sub("recombinations\\.coverage$", "recombinations", cov_base)
    candidates <- count_files[grepl(core, count_files)]
    if (length(candidates) > 0) {
      count_file <- candidates[1]
    }
  }
  
  if (is.null(count_file)) {
    message("  No matching *_count.txt file found for this coverage; skipping.")
    next
  }
  
  message("  Using count file: ", count_file)
  
  # ----- Load summarized junction/deletion data -----
  data <- read.table(count_file, sep = "", header = TRUE)
  data <- data[order(data$depth), ]
  data <- data[data$end < 29903, ]
  
  junction_depth <- sum(data$depth)
  message("  Total junction depth: ", junction_depth)
  
  # ----- Load coverage -----
  coverage <- read.table(coverage_path, sep = "", header = FALSE)
  coverage <- coverage %>%
    rename(Genome = V1, Position = V2, Depth = V3)
  
  total_depth <- sum(coverage$Depth)
  message("  Total mapped depth: ", total_depth)
  
  # ----- Forward deletions only -----
  data_forward <- data[data$start < data$end, ]
  
  data_forward$start_cov <- coverage[data_forward$start,  "Depth"]
  data_forward$end_cov   <- coverage[data_forward$end + 1, "Depth"]
  
  data_forward$MinCov <- pmin(data_forward$start_cov, data_forward$end_cov)
  data_forward$Frequency <- data_forward$depth / data_forward$MinCov
  data_forward$logFreq <- log10(data_forward$Frequency)
  data_forward$deletion_length <- data_forward$end - data_forward$start
  
  # ----- Filtering -----
  data_forward <- data_forward[
    data_forward$depth >= 5 &
      data_forward$Frequency >= 0.01 &
      data_forward$MinCov > 20 &
      data_forward$depth_positive > 2 &
      data_forward$depth_negative > 2 &
      data_forward$max_right_overhang > 34 &
      data_forward$max_left_overhang > 34,
  ]
  
  base_name <- sub("\\.coverage$", "", cov_base)
  output_table <- paste0(base_name, "_0.01_5_20_2_34.txt")
  write.table(data_forward, file = output_table, sep = "\t", row.names = FALSE)
  message("  Wrote filtered table: ", output_table)
  
  if (nrow(data_forward) == 0) {
    message("  No junctions passed filters; skipping plot.")
    next
  }
  
  plot <- ggplot(data_forward, aes(end, start, alpha = logFreq)) +
    geom_point(size = 2, color = "blue") +
    theme_linedraw(base_size = 20) +
    xlim(0, 31500) + ylim(0, 31500) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title     = element_blank(),
      legend.text      = element_text(size = 11),
      legend.position  = c(0.15, 0.63)
    ) +
    labs(x = "Stop Position (nt)", y = "Start Position (nt)") +
    scale_alpha_continuous(
      limits = c(-2, 0),
      breaks = seq(-2, 0, by = 0.5),
      labels = c("-2.0", "-1.5", "-1.0", "-0.5", "0")
    )
  
  print(plot)
  
  output_plot <- paste0(base_name, "_0.01_5_20_2_34_blue.tiff")
  ggsave(
    filename  = output_plot,
    plot      = plot,
    scale     = 1,
    width     = 4.5,
    height    = 4,
    units     = "in",
    dpi       = 1200,
    limitsize = TRUE
  )
  message("  Wrote plot: ", output_plot)
}

message("\nAll done.")
