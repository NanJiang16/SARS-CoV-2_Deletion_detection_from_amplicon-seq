#!/usr/bin/env Rscript

library(tidyverse)

# List all files ending with 'recombinations.sorted.txt' in the current working directory
input_files <- list.files(pattern = "*recombinations.sorted.txt")

# Function to process each file
process_file <- function(file) {
  # Load the data
  data <- read_tsv(file)
  
  # Summarize the data
  summarized_data <- data %>%
    group_by(chrom, start, end) %>%
    summarize(depth = n(),
              depth_positive = sum(strand == '+'),
              depth_negative = sum(strand == '-'),
              n_fragments = n_distinct(read_name),
              max_left_overhang = max(left_overhang),
              max_right_overhang = max(right_overhang),
              max_aligned_length = max(aligned_length))
  
  # Generate the output file name
  output_file <- sub("recombinations.sorted.txt$", "count.txt", file)
  
  # Write the summarized data to the output file
  write_tsv(summarized_data, file = output_file)
}

# Apply the process_file function to each input file
lapply(input_files, process_file)
