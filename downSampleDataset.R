library(Biostrings)
library(rtracklayer)
library(dplyr)
setwd('/workdir/jz963/Evolutionary_constraint/5_Poaceae/4_windows/')
fnames <- list.files(pattern = 'windows_512.chr.bed.txt')

set.seed(713)
for(i in 1:length(fnames)){
  curMat <- data.table::fread(input = fnames[i], sep = '\t', header = F, quote = '',
                              stringsAsFactors = F)
  curMat <- as.data.frame(curMat)
  prefix <- gsub(pattern = '.windows_512.chr.bed.txt', replacement = '', x = fnames[i])
  CDS <- curMat[which(curMat$V6 == 'CDS'), ]
  CDSNum <- nrow(CDS)*2
  # for upstream
  Upstream <- curMat[which(curMat$V6 == 'Upstream'), ]
  if(nrow(Upstream) >= CDSNum/2){
    Upstream <- Upstream[sample(1:nrow(Upstream), CDSNum*0.1), ]
  }
  # for downstream
  Downstream <- curMat[which(curMat$V6 == 'Downstream'), ]
  if(nrow(Downstream) >= CDSNum/2){
    Downstream <- Downstream[sample(1:nrow(Downstream), CDSNum*0.1), ]
  }
  # for Intron
  Intron <- curMat[which(curMat$V6 == 'Intron'), ]
  if(nrow(Intron) >= CDSNum/2){
    Intron <- Intron[sample(1:nrow(Intron), CDSNum*0.1), ]
  }
  # for TE
  TE <- curMat[which(curMat$V6 == 'TE'), ]
  if(nrow(TE) >= CDSNum/2){
    TE <- TE[sample(1:nrow(TE), CDSNum*0.05), ]
  }
  # for intergenic
  intergenic <- curMat[which(curMat$V6 == 'intergenic'), ]
  if(nrow(intergenic) >= CDSNum/2){
    intergenic <- intergenic[sample(1:nrow(intergenic), CDSNum*0.05), ]
  }
  
  resMat <- rbind(CDS, Upstream, Intron, Downstream, intergenic, TE)
  colnames(resMat) <- c('chr', 'start', 'end', 'width', 'strand', 'feature')
  resGR <- GenomicRanges::GRanges(resMat[,1:3])
  resGR <- resGR[which(width(resGR) == 512)]
  resMat <- resMat[which(width(resGR) == 512), ]
  write.table(x = resMat, file = paste0('../5_Datasets/', fnames[i]),
              sep = '\t', quote = F, row.names = F, col.names = T)
}

