library(data.table)
library(dplyr)
library(readr)

rm(list = ls())

output_path <- commandArgs(trailingOnly=TRUE)[1]

methylation <- fread("rapid.gene_summary.txt", sep = "\t", header = F)
methylation <- methylation[,c(4,10)]
colnames(methylation) <- c("gene", "methylation")

methylation_genome <- methylation %>% summarize(m=mean(methylation))
write_csv(as.data.frame(methylation_genome),
          output_path,
          col_names = F)

