library(tidyverse)

# command line arguments are: deletion_table depth_table
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript summarize_deletions.R read_deletion_table depth_table")
}
deletion_table_filename <- args[1]
depth_table_filename <- args[2]

# read in the deletion table
deletion_table <- read_tsv(deletion_table_filename)

# read in the depth table
depth_table <- read_tsv(depth_table_filename, col_names = c("chrom", "position", "depth"))

# modify depth table so that it is 0-based
depth_table$position <- depth_table$position - 1

# collapse deletions with the same coordinates
deletion_table %>%
    group_by(chrom, start, end) %>%
    summarize(depth = n(),
              depth_positive = sum(strand == '+'),
              depth_negative = sum(strand == '-'),
              n_fragments = n_distinct(read_name),
              max_left_overhang = max(left_overhang),
              max_right_overhang = max(right_overhang),
              max_aligned_length = max(aligned_length)) ->
    deletion_summary

# join with depth table to get coverage before and after each deletion
deletion_summary %>%
    mutate(start_cov_position = start - 1,
           end_cov_position = end) %>%
    left_join(depth_table %>% rename(start_cov = depth), 
              by = c("chrom" = "chrom", "start_cov_position" = "position")) %>%
    left_join(depth_table %>% rename(end_cov = depth),
              by = c("chrom" = "chrom", "end_cov_position" = "position")) %>%
    mutate(min_cov = pmin(start_cov, end_cov),
           frequency = depth / min_cov) %>%
    select(-start_cov_position, -end_cov_position) ->
    deletion_summary_with_depth

# write out the summary to stdout
write_tsv(deletion_summary_with_depth, stdout())
