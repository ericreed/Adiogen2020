require(edgeR)
require(pheatmap)
require(openxlsx)
require(parallel)
require(Biobase)
require(limma)
source("parameters.R")

# Set script tag, so we know what script outputs came from
sTag <- "032a"

# Main directories

# Files to write out
outGSfile <- "LIMMA_DGE_expA"

# Begin

# Read in expression set
eSet <- readRDS(eSetCombatFile011a)

# Run LIMMA
outDA <- runLIMMA(eSet, control = "Vehicle Control", outDir = gsigDir, sTag = sTag, outGSfile, cores = 4)
writeXLSX(outDA, LIMMAoutFile032a)

# Subset for signficant genes (FDR < 0.1)
outSig <- lapply(outDA, function(x) x[x$adj.P.Val < 0.1,])
writeXLSX(outSig, LIMMAoutFileSig032a)

# Create dataMatrix of test statistics
Dlist <- lapply(outDA, function(x){
  x$z <- qnorm(x$P.Value/2, lower.tail = FALSE) * sign(x$logFC)
  rownames(x) <- toupper(x$gene)
  return(x)
})
names(Dlist) <- names(outDA)

# Create matrix of test statistics
genes <- rownames(Dlist[[1]])
Dfram <- as.data.frame(do.call(cbind, lapply(Dlist, function(x) x[genes, "z"] ) )); rownames(Dfram) <- genes

# Fix chemical abbs
update <- read.csv(updateFile, stringsAsFactors = F)
update$Abbreviation <- toupper(update$Abbreviation)
rownames(update) <- update$Abbreviation
colnames(Dfram) <- update[colnames(Dfram),"Update"]

# Write matrix
saveRDS(Dfram, dataMatrixFile032a)

# Empty Environment
rm(list = ls())
