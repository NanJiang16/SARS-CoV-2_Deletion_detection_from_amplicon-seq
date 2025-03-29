#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

# -----------------------------
# STEP 1: Summarize recombination data
# -----------------------------
recomb_files <- list.files(pattern = "*recombinations.sorted.txt")

process_file <- function(file) {
  data <- read_tsv(file)
  
  summarized_data <- data %>%
    group_by(chrom, start, end) %>%
    summarize(depth = n(),
              depth_positive = sum(strand == '+'),
              depth_negative = sum(strand == '-'),
              n_fragments = n_distinct(read_name),
              max_left_overhang = max(left_overhang),
              max_right_overhang = max(right_overhang),
              max_aligned_length = max(aligned_length)) %>%
    ungroup()
  
  output_file <- sub("recombinations.sorted.txt$", "count.txt", file)
  write_tsv(summarized_data, file = output_file)
}

lapply(recomb_files, process_file)

# -----------------------------
# STEP 2: Frequency calculation and filtering
# -----------------------------
coverage_files <- list.files(path = ".", pattern = "*_recombinations.coverage", full.names = TRUE)
bed_files <- list.files(path = ".", pattern = "*_count.txt", full.names = TRUE)

for(i in seq_along(coverage_files)){
  coverage_file <- coverage_files[i]
  base_name <- gsub("_recombinations\\.coverage$", "", basename(coverage_file))
  bed_pattern <- paste0(base_name, "_count\\.txt")
  bed_file <- bed_files[grepl(pattern = bed_pattern, bed_files)]
  
  if(length(bed_file) == 0){
    print(paste("No count file found for:", base_name))
    next
  }
  
  data <- read.table(bed_file, sep = "", header = TRUE)
  data <- data[order(data$depth), ]
  data <- data[which(data$end < 29903), ]
  junction_depth <- sum(data$depth)
  print(paste0("The total junction depth is: ", junction_depth))
  
  coverage <- read.table(coverage_file, sep = "", header = FALSE)
  coverage <- coverage %>% rename(Genome = V1, Position = V2, Depth = V3)
  total_depth <- sum(coverage$Depth)
  print(paste0("The total mapped depth is: ", total_depth))
  
  data_forward <- data[which(data$start < data$end), ]
  data_forward$start_cov = coverage[unlist(data_forward[2]), 3]
  data_forward$end_cov = coverage[unlist(data_forward[3] + 1), 3]
  data_forward <- transform(data_forward, MinCov = pmin(start_cov, end_cov))
  data_forward$Frequency = data_forward$depth / data_forward$MinCov
  data_forward$logFreq = log10(data_forward$Frequency)
  data_forward$deletion_length = data_forward$end - data_forward$start
  
  data_forward <- data_forward[
    data_forward$depth >= 5 &
      data_forward$Frequency >= 0.01 &
      data_forward$MinCov > 20 &
      data_forward$depth_positive > 2 &
      data_forward$depth_negative > 2 &
      data_forward$max_right_overhang > 34 &
      data_forward$max_left_overhang > 34,
  ]
  
  write.table(data_forward, file = paste0(base_name, "_0.01_5_20_2_34.txt"), sep = "\t", row.names = FALSE)
  
  plot <- ggplot(data_forward, aes(end, start, alpha = logFreq)) + 
    geom_point(size = 2, color = "blue") + 
    theme_linedraw(base_size = 20) + 
    ylim(0, 31500) + xlim(0, 31500) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          legend.position = c(0.15, 0.63)) + 
    labs(x = "Stop Position (nt)", y = "Start Position (nt)") +
    scale_alpha_continuous(limits = c(-2, 0), 
                           breaks = seq(-2, 0, by = 0.5), 
                           labels = c("-2.0", "-1.5", "-1.0", "-0.5", "0"))
  
  print(plot)
  ggsave(filename = paste0(base_name, "_0.01_5_20_2_34_blue.tiff"), plot = plot,
         scale = 1, width = 4.5, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
}

