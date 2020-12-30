require(edgeR)
require(pheatmap)
require(Biobase)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
source("parameters.R")

# Set script tag, so we know what script outputs came from
sTag <- "003"

# Set plate ID
pID <- "plate4"
pNo <- 4


## Subdirectories
platDir <- file.path(datDir, "bulk", pID)

# Files to read in
platFile <- file.path(platDir, "SSF-2156.HHKY3BGX2.unq.refseq.umi.dat")

# Files to write out

## Begin

## Read in data files

### Sample Information
info <- read.csv(SampFile, stringsAsFactors = FALSE)
info <- info[info$Plate == pNo,] # Subset for plate 1
row.names(info) <- info$Position # We will match the data based on position on the plate

### Expression data
plat <- read.table(platFile, sep = "\t", header = T, row.names = 1)
colnames(plat) <- sub("[[:alnum:]]+_", "", colnames(plat)) # Get position
info <- info[colnames(plat),]

### Keep RefSeq genes only
plat <- plat[!grepl("Rik", row.names(plat)),]

### Create expression set
eSet <- to.eSet(plat, info, data.frame(RefSeq = row.names(plat), row.names = row.names(plat)))
saveRDS(eSet, eSetFile003)
### Create a log2CPM transformed eSet
eSetCPM <- eSet
exprs(eSetCPM) <- cpm(eSetCPM, log = TRUE)

### Boxplots Per Sample (CPM)
creatBoxPlots <- function(eset, title){
  e <- exprs(eset)
  i <- data.frame(Chemical = pData(eset)$Chemical, Experiment =  pData(eset)$Experiment, Sample = row.names(pData(eset)))
  
  eG <- melt(e)
  colnames(eG) <- c("Gene", "Sample", "Expression")
  eG <- merge(eG, i)
  eG <- eG[order(eG$Chemical),]
  
  p <- ggplot(eG, aes(Chemical, Expression, fill = Sample)) +
    geom_boxplot() +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title)
  return(p)
}

png(AboxplotFile003, width = 1200, height = 600)
creatBoxPlots(eSetCPM[,pData(eSetCPM)$Experiment == "A"], title =  "Experiment A (log2CPM)")
dev.off()

png(BboxplotFile003, width = 1200, height = 600)
creatBoxPlots(eSetCPM[,pData(eSetCPM)$Experiment == "B"], title =  "Experiment B (log2CPM)")
dev.off()

### Barplots Per Sample (Sum of Counts)
creatBarPlots <- function(eset, title){
  e <- exprs(eset)
  i <- data.frame(Chemical = pData(eset)$Chemical, Experiment =  pData(eset)$Experiment, Sample = row.names(pData(eset)))
  
  eG <- data.frame(Sample = colnames(e), log2Counts = log(colSums(e), 2))
  eG <- merge(eG, i)
  eG <- eG[order(eG$Chemical),]
  
  p <- ggplot(eG, aes(x = Chemical, y = log2Counts, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(legend.position = "None",
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title)
  return(p)
}

png(AbarplotFile003, width = 1200, height = 600)
creatBarPlots(eSet[,pData(eSet)$Experiment == "A"], title =  "Experiment A (Total Reads)")
dev.off()

png(BbarplotFile003, width = 1200, height = 600)
creatBarPlots(eSet[,pData(eSet)$Experiment == "B"], title =  "Experiment B (Total Reads)")
dev.off()

### Check relationship between expression and mad score
madVec <- apply(exprs(eSetCPM), 1, mad)
medVec <- apply(exprs(eSetCPM), 1, median)
madVmed <- data.frame(MADscore = madVec, MedExp = medVec)
png(medVmadplotFile003, 400, 400)
ggplot(madVmed, aes(MedExp, MADscore)) + geom_point()
dev.off()

### Cluster top 3K mad genes
madVec <- sort(madVec, decreasing = TRUE)
eSetMAD <- eSetCPM[names(madVec)[1:3000],]

png(clustMADplotFile003, width = 1000, height = 1000)
MADclust <- clusEset(eSetMAD)
dev.off()

### Cluster top 3K Expressed genes
medVec <- sort(medVec, decreasing = TRUE)
eSetMED <- eSetCPM[names(medVec)[1:3000],]

png(clustExpplotFile003, width = 1000, height = 1000)
MEDclust <- clusEset(eSetMED)
dev.off()

# Save
eSetRem <- eSet
saveRDS(eSetRem, eSetRemFile003)

# Empty Environment
rm(list = ls())