library(Biostrings)
library(dplyr)
library(ggplot2)

# Function to process a single file
process_file <- function(file) {
  SARS2_TRS = "ACGAAC"
  SARS2_fasta <- readDNAStringSet("./NC_045512.2.fasta", "fasta")
  TRS_positions <- vmatchPattern(SARS2_TRS, SARS2_fasta)
  TRS_start <- startIndex(TRS_positions)
  TRS_stop <- endIndex(TRS_positions)
  TRS_start_window <- sapply(TRS_start, function(x) x - 31) #set window for TRS start positions
  TRS_stop_window <- sapply(TRS_stop, function(x) x + 31) #set window for TRS stop positions
  colnames(TRS_start_window) <- c("Start")
  colnames(TRS_stop_window) <- c("Stop")
  TRS_matrix <- merge(TRS_start_window, TRS_stop_window, by="row.names")
  TRS_matrix <- as.matrix(TRS_matrix[-1])
  if(nrow(TRS_matrix) > 9) {
    stop("Too many TRS positions detected in the genome. Virus may be evolving new TRS locations. Script aborting.")
  }
  if(nrow(TRS_matrix) < 9) {
    stop("Not enough TRS sequences detected in the genome. Virus may be using a mutated version of the consensus sequence. Script aborting.")
  }
  df <- read.table(file, header = TRUE)
  df_TRSL <- filter(df, between(df$start, as.numeric(TRS_matrix[1,1]), as.numeric(TRS_matrix[1,2])))
  df_sgmRNA2 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[2,1]), as.numeric(TRS_matrix[2,2])))
  df_sgmRNA3 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[3,1]), as.numeric(TRS_matrix[3,2])))
  df_sgmRNA4 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[4,1]), as.numeric(TRS_matrix[4,2])))
  df_sgmRNA5 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[5,1]), as.numeric(TRS_matrix[5,2])))
  df_sgmRNA6 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[6,1]), as.numeric(TRS_matrix[6,2])))
  df_sgmRNA7 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[7,1]), as.numeric(TRS_matrix[7,2])))
  df_sgmRNA8 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[8,1]), as.numeric(TRS_matrix[8,2])))
  df_sgmRNA9 <- filter(df_TRSL, between(df_TRSL$end, as.numeric(TRS_matrix[9,1]), as.numeric(TRS_matrix[9,2])))
  #Create null dataframes if no observations
  if(nrow(df_sgmRNA2) == 0) {
    df_sgmRNA2 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA3) == 0) {
    df_sgmRNA3 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA4) == 0) {
    df_sgmRNA4 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA5) == 0) {
    df_sgmRNA5 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA6) == 0) {
    df_sgmRNA6 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA7) == 0) {
    df_sgmRNA7 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA8) == 0) {
    df_sgmRNA8 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  if(nrow(df_sgmRNA9) == 0) {
    df_sgmRNA9 <- data.frame("chrom" = "MN908947.3", "start" = 0, "end" = 0, "depth" = 0, "depth_positive" = 0, "depth_negative" = 0, "n_fragments" = 0, "max_left_overhang" = 0, "max_right_overhang" = 0, "max_aligned_length"=0, "start_cov"=0, "end_cov"=0, "MinCov"=0, "Frequency"=0, "logFreq"=0, "deletion_length"=0)
  }
  #Add column identifying sgmRNA species
  if(nrow(df_sgmRNA2) > 0) {
    df_sgmRNA2$Type <- "sgmRNA2"
  }
  
  if(nrow(df_sgmRNA3) > 0) {
    df_sgmRNA3$Type <- "sgmRNA3"
  }
  if(nrow(df_sgmRNA4) > 0) {
    df_sgmRNA4$Type <- "sgmRNA4"
  }
  
  if(nrow(df_sgmRNA5) > 0) {
    df_sgmRNA5$Type <- "sgmRNA5"
  }
  if(nrow(df_sgmRNA6) > 0) {
    df_sgmRNA6$Type <- "sgmRNA6"
  }
  
  if(nrow(df_sgmRNA7) > 0) {
    df_sgmRNA7$Type <- "sgmRNA7"
  }
  if(nrow(df_sgmRNA8) > 0) {
    df_sgmRNA8$Type <- "sgmRNA8"
  }
  
  if(nrow(df_sgmRNA9) > 0) {
    df_sgmRNA9$Type <- "sgmRNA9"
  }
  #Slice canonical sgmRNA species
  slice <- dplyr::slice
  sgmRNA2_canonical <- df_sgmRNA2 %>% arrange(desc(df_sgmRNA2$depth)) %>% slice(1)
  sgmRNA3_canonical <- df_sgmRNA3 %>% arrange(desc(df_sgmRNA3$depth)) %>% slice(1)
  sgmRNA4_canonical <- df_sgmRNA4 %>% arrange(desc(df_sgmRNA4$depth)) %>% slice(1)
  sgmRNA5_canonical <- df_sgmRNA5 %>% arrange(desc(df_sgmRNA5$depth)) %>% slice(1)
  sgmRNA6_canonical <- df_sgmRNA6 %>% arrange(desc(df_sgmRNA6$depth)) %>% slice(1)
  sgmRNA7_canonical <- df_sgmRNA7 %>% arrange(desc(df_sgmRNA7$depth)) %>% slice(1)
  sgmRNA8_canonical <- df_sgmRNA8 %>% arrange(desc(df_sgmRNA8$depth)) %>% slice(1)
  sgmRNA9_canonical <- df_sgmRNA9 %>% arrange(desc(df_sgmRNA9$depth)) %>% slice(1)
  #Create concatenated dataframe of canonical sgmRNAs
  df_canonical <- rbind(sgmRNA2_canonical, sgmRNA3_canonical, sgmRNA4_canonical, sgmRNA5_canonical, sgmRNA6_canonical, sgmRNA7_canonical, sgmRNA8_canonical, sgmRNA9_canonical)
  df_canonical$Total <- sum(df_canonical$depth)
  df_sgmRNA <- rbind(df_sgmRNA2, df_sgmRNA3, df_sgmRNA4, df_sgmRNA5, df_sgmRNA6, df_sgmRNA7, df_sgmRNA8, df_sgmRNA9)
  df_sgmRNA$Total <- sum(df_sgmRNA$depth)
  #Print list of alternative sgmRNAs
  df_alternative <- anti_join(df_sgmRNA, df_canonical, by = c("start", "end"))
  if(nrow(df_alternative) > 0) {
    df_alternative$Total <- sum(df_alternative$depth)
  }
  df_alt_summary <- df_alternative %>% group_by(Type) %>% summarise(Sum = sum(depth))
  #Slice DVGs and turn format into BED
  df_DVG <- anti_join(df, df_sgmRNA, by = c("start", "end"))
  #save tables and scatter plots
  write.table(df_canonical, file = paste0(gsub(".txt", "", file), "_canonical_sgmRNAs.txt"), sep = "\t", row.names = FALSE)
  write.table(df_alternative, file = paste0(gsub(".txt", "", file), "_alternative_sgmRNAs.txt"), sep = "\t", row.names = FALSE)
  write.table(df_sgmRNA, file = paste0(gsub(".txt", "", file), "_sgmRNAs.txt"), sep = "\t", row.names = FALSE)
  write.table(df_alt_summary, file = paste0(gsub(".txt", "", file), "_alt_sgmRNA_summary.txt"), sep = "\t", row.names = FALSE)
  write.table(df_DVG, file = paste0(gsub(".txt", "", file), "_DVGs.bed.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  data_graph1 <- ggplot(df_DVG, aes(end, start, alpha = logFreq)) + 
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
          scale_alpha_continuous(limits = c(-5, 0), 
                                 breaks = seq(-5, 0, by = 1), 
                                 labels = c("-5.0", "-4.0", "-3.0", "-2.0", "-1.0", "0"))
  )
  
}
# List all .txt files in the current directory
files <- list.files(pattern = "*_65_20_1_75_0.01_5_20_2_34.txt")

# Apply the function to each file
lapply(files, process_file)

