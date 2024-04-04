library(rtracklayer)
library(snowfall)
library(GenomicFeatures)
library(Biostrings)
library(argparse)
library(parallel)
options(scipen = 20)

# print usage of the script
printUsage <- function(){
  cat("Usage: Rscript genomic_windows_class.R <input BED file> <GFF file> <TE file>\n")
}

# get overlap with each feature
getPercentage <- function(queryGR, subjectGR){
  hits <- findOverlaps(query = queryGR, subject = subjectGR)
  overlaps <- pintersect(queryGR[queryHits(hits)], subjectGR[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(queryGR[queryHits(hits)])
  resDF <- data.frame(ID = 1:length(queryGR))
  resDF$percentage <- 0
  resDF$percentage[queryHits(hits)] <- percentOverlap
  resDF
}
args=commandArgs(T)

fname <- args[1]
gffName <- args[2]
TEName <- args[3]
threads <- as.numeric(args[4])
output <- args[5]
BEDs <- import(fname)

txdb <- makeTxDbFromGFF(file = gffName)
TE <- import(TEName)


CDSbyTx <- unique(unlist(cdsBy(x = txdb, by = c('tx', 'gene'), use.names = TRUE)))
CDSbyTx$transcript <- names(CDSbyTx)


CDSbyTx <- cdsBy(x = txdb, by = c('tx', 'gene'), use.names = TRUE)
result <- mclapply(names(CDSbyTx), function(i) {
  range_i <- as.data.frame(range(CDSbyTx[i]))
  range_i$transcript <- i
  range_i
}, mc.cores = threads)
result <- do.call(what = rbind, args = result)

CDSbyTx <- GenomicRanges::GRanges(unique(result))
rm(result)
gc()

upstream <- promoters(x = CDSbyTx, upstream = 1500, downstream = 0)
downstream <- flank(x = CDSbyTx, start = FALSE, width = 1500)

CDS <- cds(x = txdb, columns=c("cds_id", "tx_name"))
introns <- intronicParts(txdb = txdb)

cdsPerc <- getPercentage(queryGR = BEDs, subjectGR = CDS)
intronsPerc <- getPercentage(queryGR = BEDs, subjectGR = introns)
upstreamPerc <- getPercentage(queryGR = BEDs, subjectGR = upstream)
downstreamPerc <- getPercentage(queryGR = BEDs, subjectGR = downstream)
TEPerc <- getPercentage(queryGR = BEDs, subjectGR = TE)

resDF <- data.frame(cdsPerc = cdsPerc$percentage,
                    intronsPerc = intronsPerc$percentage,
                    upstreamPerc = upstreamPerc$percentage,
                    downstreamPerc = downstreamPerc$percentage,
                    TEPerc = TEPerc$percentage)

type <- c('CDS', 'Intron', "Upstream", 'Downstream', 'TE')
resType <- apply(resDF, 1, FUN = function(x){
  if(all(as.numeric(x) == 0)){
    res <- 'intergenic'
  }else{
    res <- type[which.max(as.numeric(x))]
  }
  return(res)
})

BEDs <- data.frame(BEDs, stringsAsFactors = F)
BEDs$type <- resType
write.table(x = BEDs, file = output, sep = '\t', quote = F,
            row.names = F, col.names = F)


