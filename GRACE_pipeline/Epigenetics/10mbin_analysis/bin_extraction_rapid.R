library(data.table)
library(tidyverse)

rm(list = ls())

output_path <- commandArgs(trailingOnly=TRUE)[1]

methylation_summary <- fread("rapid.10mbin_summary.txt", 
                             sep = "\t", 
                             header = F)
methylation_summary <- methylation_summary[,c(4,8)]
colnames(methylation_summary) <- c("bin", "methylation")

methylation_bin <- methylation_summary %>% group_by(bin) %>% summarize(methy=mean(methylation))

all_bin <- read.csv("Epigenetics/10mbin_analysis/all_bin.csv", 
                    stringsAsFactors = FALSE, 
                    header = F)

rapid_bin <- as.data.frame(t(methylation_bin))

rapid_bin_names <- as.character(rapid_bin[1, ])
rapid_bin_values <- as.numeric(rapid_bin[2, ])

all_bin_names <- as.character(all_bin[1, ])

for (i in 1:length(all_bin_names)) {
  match_idx <- which(rapid_bin_names == all_bin_names[i])
  
  if (length(match_idx) > 0) {
    all_bin[2, i] <- rapid_bin_values[match_idx]
  } else {
    all_bin[2, i] <- NA
  }
}

write.table(all_bin, 
            file = output_path, 
            sep = ",",
            row.names = FALSE,
            col.names = FALSE,
            na = "",
            quote = FALSE)
