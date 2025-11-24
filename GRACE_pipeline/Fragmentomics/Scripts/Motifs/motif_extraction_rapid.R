library(data.table)
library(tidyr)
library(readr)

rm(list = ls())

sample_info_path <- commandArgs(trailingOnly=TRUE)[1]
motifname_path <- commandArgs(trailingOnly=TRUE)[2]
output_path <- commandArgs(trailingOnly=TRUE)[3]

get_matrix <- function(samples_info, motifs){
  c=0
  while (c < length(rownames(samples_info))){
    c=c+1
    path <- samples_info[c,1]
    sname <- samples_info[c,2]
    file_data <- as.data.frame(readRDS(file.path(path,paste(sname,".motif.R",sep=""))))
    file_data <- subset(file_data, file_data$Var1 %in% motifs$End.motif)
    colnames(file_data)[2] <- sname
    file_data[,2] <- as.numeric(file_data[,2]) 
    file_data <- cbind(sname, file_data)
    file_data[,4] <- file_data[,3]/sum(file_data[,3])*100
    colnames(file_data) <- c("samplename", "motif","count", "normalizedToCounts")
    
    if (c == 1){
      total <- file_data
    } else {
      total <- rbind(total, file_data)
    } 
  }
  return(total)
}

#### load sample annotation file
samples_info <- read.table(sample_info_path,
                           sep="\t",
                           header=F,
                           stringsAsFactors = F)

#### select motif order
motifNames <- read.delim(motifname_path,
                         header = T,
                         sep = "\t")

#### create matrix from sample
df <- get_matrix(samples_info,
                 motifNames)

#### create matrix of "normalizedToCounts" values
groupby = "normalizedToCounts"
dt<-setDT(df,
          key = "samplename")  
perSample = data.table(motif=unique(dt$motif),
                       key="motif")
for (s in unique(df$samplename)) {
  subTab = dt[as.character(s)]
  groupby=paste0("^",groupby,"$")
  ind=grep(groupby,colnames(dt),value = T)
  perSample[as.character(subTab$motif),as.character(s)] = subTab[, ..ind]
}
perSample<-as.data.frame(perSample)

#### save significant motifs
motif_qtop10 <- as.vector(c("CAAG","CACC",	"CTCA",	"CTCC",	"CTCT",	"CTTG",	"GCAA",	"GCAT",	"TACG",	"TTGA"))

motif_matrix <- perSample[perSample$motif %in% motif_qtop10,]
write_csv(as.data.frame(t(motif_matrix)),
          output_path,
          col_names = F)



