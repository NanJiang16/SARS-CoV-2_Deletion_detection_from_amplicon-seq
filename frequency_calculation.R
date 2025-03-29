#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

# get list of coverage files
coverage_files <- list.files(path = ".", pattern = "*_recombinations.coverage", full.names = TRUE)

# get list of bed files
bed_files <- list.files(path = ".", pattern = "*_count.txt", full.names = TRUE)

# loop over each coverage file
for(i in seq_along(coverage_files)){
  coverage_file <- coverage_files[i]
  
  # Extract the common base name by removing the '_concatenate_recombinations.coverage' suffix
  base_name_pattern <- gsub("_recombinations\\.coverage$", "", basename(coverage_file))
  
  # Create a pattern to find the matching bed file by replacing the coverage suffix with the bed suffix
  bed_pattern <- paste0(base_name_pattern, "_count\\.txt")
  
  # Find the paired bed file using the pattern
  bed_file <- bed_files[grepl(pattern = bed_pattern, bed_files)]
  
  # Check if corresponding bed file exists
  if(length(bed_file) == 0){
    print(paste("No bed file found for base name:", base_name_pattern))
    next
  }
  
  # now replace all instances of your hardcode file names with variables
  data <- read.table(bed_file, sep = "", header = TRUE)
  data <- data[,c(1,2,3,4,5,6,7,8,9,10)]
  data <- data[order(data$depth), ]
  data <- data[which(data$end < 29903), ]
  junction_depth <- sum(data$depth)
  print(paste0("The total junction depth is: ", junction_depth))
  
  coverage <- read.table(coverage_file, sep="", header = FALSE)
  coverage <- coverage %>% rename(Genome = V1, Position = V2, Depth = V3)
  
  total_depth <- sum(coverage$Depth)
  print(paste0("The total mapped depth is: ", total_depth))
  
  data_forward <- data[which(data$start < data$end), ]
  data_forward$start_cov = coverage[unlist(data_forward[2]),3]
  data_forward$end_cov = coverage[unlist(data_forward[3]+1),3]
  data_forward <- transform(data_forward, MinCov = pmin(start_cov, end_cov))
  data_forward$Frequency = data_forward$depth / data_forward$MinCov
  data_forward$logFreq = log10(data_forward$Frequency)
  data_forward$deletion_length = data_forward$end - data_forward$start
  data_forward <- data_forward[data_forward$depth >= 5 & data_forward$Frequency >= 0.01 & data_forward$MinCov > 20 & data_forward$depth_positive > 2 & data_forward$depth_negative > 2 & data_forward$max_right_overhang > 34 & data_forward$max_left_overhang > 34, ]
  
  
  write.table(data_forward, file = paste0(base_name_pattern, "_0.01_5_20_2_34.txt"), sep = "\t", row.names = FALSE)
  data_graph1 <- ggplot(data_forward, aes(end, start, alpha = logFreq)) + 
    geom_point(size = 2, color = "blue") + 
    theme_linedraw(base_size = 20) + 
    ylim(0, 31500) + xlim(0, 31500) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.position = c(0.15, 0.63)) + 
    labs(x = "Stop Position (nt)", y = "Start Position (nt)")
  
  print(data_graph1 + 
          scale_alpha_continuous(limits = c(-2, 0), 
                                 breaks = seq(-2, 0, by = 0.5), 
                                 labels = c("-2.0", "-1.5", "-1.0", "-0.5", "0"))
  )
  ggsave(filename = paste0(base_name_pattern, "_0.01_5_20_2_34_blue.tiff"), plot = last_plot(), scale = 1, width = 4.5, height = 4, units = "in", dpi = 1200, limitsize = TRUE)
  
}