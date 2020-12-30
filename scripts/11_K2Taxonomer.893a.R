require(pheatmap)
require(openxlsx)
require(dplyr)
require(tidyr)
require(limma)
require(dendextend)
require(ggdendro)
require(ggplot2)
require(Biobase)
require(RColorBrewer)
require(cba)
source("parameters.R")

# Set script tag, so we know what script outputs came from
sTag <- "893"

# Main directories

## Subdirectories

# Read in eSet
eSet <- readRDS(eSetCombatFile011a)
rownames(eSet) <- toupper(rownames(eSet))

# Files to Read in
gSigListFiles <- list.files(gsigDir, pattern = "032a.rds")

# Files to write out

# Read in Differential Results
Dlist <- lapply(gSigListFiles, function(x){
  x <- readRDS(file.path(gsigDir, x))@Signature
  x$Z <- qnorm(x$P.Value/2, lower.tail = FALSE)*sign(x$Score)
  x$FDR <- p.adjust(x$P.Value, method = "BH")
  x$ID <- toupper(x$ID)
  rownames(x) <- x$ID
  return(x)
})
names(Dlist) <- gsub("LIMMA_DGE_expA_|_Vehicle_032a.rds", "", gSigListFiles)

# Only Use True or Predicted PPARg modifiers
PgMod <- unique(eSet$Chemical[eSet$PPARg_Mod == 'Yes'])
PgMod <- PgMod[!is.na(PgMod)]
PgMod <- c(PgMod, "TONALID", "QUINOXYFEN", "ALLETHRIN", "FENTHION", "NAIVE")

Dlist <- Dlist[PgMod]

# Create matrix of test statistics
genes <- rownames(Dlist[[1]])
Dfram <- as.data.frame(do.call(cbind, lapply(Dlist, function(x) x[genes, "Z"] ) )); rownames(Dfram) <- genes

# Fix chemical abbs
update <- read.csv(updateFile, stringsAsFactors = F)
update$Abbreviation <- toupper(update$Abbreviation)
rownames(update) <- update$Abbreviation
colnames(Dfram) <- update[colnames(Dfram),"Update"]

# Add vehicle
Dfram$VEHICLE <- 0

# Create list of splits comparing agg methods

## Ward's Method
finUnListWard <- runK2Clust(Dfram, nGenes = nrow(Dfram)/10, agg = "ward.D2")
splitNames <- paste0(LETTERS[1:length(finUnListWard)], "(", unlist(lapply(finUnListWard, function(x) format(round(x$bootP, 2), nsmall = 2))), ")")
dClustWard <- mat2dendro(finUnListWard, splitNames = splitNames)

# Which to use
finUnList <- finUnListWard

# Create dendrogram
dClust <- mat2dendro(finUnList, splitNames = splitNames)
plot(dClust)

# Compare to regular clustering
RegClust <- runHClust(Dfram, nGenes = nrow(Dfram)/10)
png(HclustPlotOut, width = 1000, height = 1100)
plot_horiz.dendrogram(RegClust)
dev.off()

# Plot dendrogram
# Add Labels (of segments)
png(dendroPlotOut, width = 900, height = 1000)
dendroP <- plotDendro(dClust, splitNames = splitNames)
dendroP
dev.off()

# Create dendro as hierarchichal network
vNetOut <- generateVisNetwork(finUnList)

# Add DGE results to list of splits
# Fix chemical IDs
eSet$Chemical <- update[eSet$Chemical,"Update"]
finUnList <- runDGE_mods(eSet, finUnList, Dfram, nGenes = nrow(Dfram)/10)

# Fix FDRs
finUnList <- fixFDRs(finUnList)

# Get gene lists for each split
names(finUnList)  <- LETTERS[1:length(finUnList)]
gLists <- lapply(finUnList, function(x){
  x <- x$res
  x <- x[x$adj.P.Val < 0.1 & abs(x$logFC) > 1,]
  x$diff <- apply(abs(x[, c("mean_z1", "mean_z2")]), 1, which.max)
  list(up1 = rownames(x)[x$diff == 1 & x$mean_z1 > 0],
       down1 = rownames(x)[x$diff == 1 & x$mean_z1 < 0],
       up2 = rownames(x)[x$diff == 2 & x$mean_z2 > 0],
       down2 = rownames(x)[x$diff == 2 & x$mean_z2 < 0])
})


# Get set of top gene markers for plotting purposes
terminalsplits <- c("G", "J", "I", "K", "L", "N", "M", "E", "H", "D")
topGenes <- do.call(rbind, lapply(seq(length(terminalsplits)), function(x){
  let <- terminalsplits[x]
  xSub <- finUnList[[let]]$res[finUnList[[let]]$res$adj.P.Val < 0.1 & abs(finUnList[[let]]$res$logFC) > 1, ]
  xSub$cluster <- let
  xSub$clustOrd <- x
  xSub$group <- c(1, 2)[as.numeric(abs(xSub$mean_z2) > abs(xSub$mean_z1)) + 1]
  return(xSub)
}))

# Keep only splits that lead to terminal ndoes
termList <- list(c(2, 1), c(2, 1), c(2), c(2), c(2), c(2, 1), c(1), c(2, 1), c(2, 1), c(1))
names(termList) <- terminalsplits
topGenes <- do.call(rbind, lapply(terminalsplits, function(let) {
  topGenes[topGenes$cluster == let & topGenes$group %in% termList[[let]],]
}))

# Order by p-value
topGenes <- topGenes[order(topGenes$P.Value),]

# Get unique gene lists
topGenesUnique <- topGenes[!duplicated(topGenes$RefSeq),]

# Subset Dfram
# Cluster each gene set
DframSort <- do.call(rbind, lapply(terminalsplits, function(clustID){
  # Subset
  topSub <- topGenesUnique[topGenesUnique$cluster == clustID,]

  DframClust <- do.call(rbind, lapply(rev(sort(unique(topSub$group))), function(g){
    genes <- toupper(as.character(topSub$RefSeq[topSub$group == g]))
    DframSub <- Dfram[genes,]
    if(nrow(DframSub) > 1){
      gDist <- as.dist(1-cor(t(Dfram[genes,])))
      gClust <- hcopt(gDist, method = "ward.D")
      DframSub <- DframSub[gClust$order,]
    }
  }))
  return(DframClust)
}))
DframSig <- DframSort[,rev(get_leaves_attr(dClust, "label"))]

# Add row annotation
annotRows <- do.call(cbind, lapply(terminalsplits, function(clustID){
  # Subset
  topSub <- topGenes[topGenes$cluster == clustID,]
  annotClust <- as.data.frame(do.call(cbind, lapply(rev(sort(unique(topSub$group))), function(g){
    c("No", "Yes")[as.numeric(rownames(DframSig) %in% toupper(as.character(topSub$RefSeq[topSub$group == g]))) + 1]
  })))
  colnames(annotClust) <- paste0(clustID, rev(sort(unique(topSub$group))))
  return(annotClust)
}))
rownames(annotRows) <- rownames(DframSig)
annotRows <- annotRows[,ncol(annotRows):1]

# Add colors to annotation
colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff')

# Add column annotations
termClusts <- unlist(lapply(terminalsplits, function(clustID) {
  chems <- finUnList[[clustID]]$chems[termList[[clustID]]]
}), recursive = F)
annotCols <- data.frame(terminal = unlist(sapply(seq(length(termClusts)), function(i) rep(colors[i], length(termClusts[[i]])))))
rownames(annotCols) <- unlist(termClusts)

## Rows
annotColors <- lapply(colors, function(x) c(No = "white", Yes = x))[rev(seq(ncol(annotRows)))]
names(annotColors) <- colnames(annotRows)

## Cols
annotColors$terminal <- rev(colors[seq(ncol(annotRows))])
names(annotColors$terminal) <- annotColors$terminal

# Generate heatmap
png(file.path(plotDir, "Tax_heatmap_893.png"), width = 1030, height = 1500)
pheatmap(DframSig,
  color = rev(brewer.pal(7, "RdBu")),
  breaks = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 20,
  show_rownames = F,
  legend = FALSE,
  annotation_legend = FALSE,
  annotation_col = annotCols,
  annotation_row = annotRows,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_colors = annotColors)
dev.off()

# Run hyperenrichment
gUnLists <- unlist(gLists, recursive = F)

# Write to file
gLists <- createWorkbook()
lapply(names(gUnLists), function(x){
  gL <- gUnLists[[x]]
  addWorksheet(gLists, sheetName = x)
  writeData(gLists, sheet = x, gL)
})
saveWorkbook(gLists, file = GeneListTableOut, overwrite = T)

# Get annotations
REAC <- getgeneLists(goFile)
REAC <- lapply(REAC, toupper)
PS273 <- getgeneLists(gmtOut_GSE22033)
REAC <- unlist(list(REAC, PS273), recursive = F)

# Find complete set of genes before filtering
eSetRaw <- readRDS(eSetAFile010a)
rownames(eSetRaw) <- toupper(rownames(eSetRaw))
REAC <- lapply(REAC, function(x) x[x %in% rownames(eSetRaw)])
REAC <- REAC[unlist(lapply(REAC, function(x) length(x) >= 10 & length(x) <= 1000))]

# Run HE
HyperEnr <- as.data.frame(hyperEnrichment(gUnLists, REAC, min.drawsize = 1, ntotal = nrow(eSetRaw)))
for(i in 1:ncol(HyperEnr)) HyperEnr[,i] <- as.character(HyperEnr[,i])
for(i in 2:7) HyperEnr[,i] <- as.numeric(HyperEnr[,i])
HyperEnr <- HyperEnr[order(HyperEnr$pval),]

# Get sig results
HEsig <- HyperEnr[HyperEnr$fdr < 0.1,]
sets <- sort(unique(HEsig$set))
setSig <- lapply(sets, function(x) HEsig$category[HEsig$set == x])
names(setSig) = sets


# Add set information for each gene
gene2Pathway <- sapply(unique(unlist(REAC)), function(x) paste(names(REAC)[unlist(lapply(REAC, function(y) x %in% y))], collapse = "; "))

# Create workbook
DGEwb <- createWorkbook()
print(dendroP)
addWorksheet(DGEwb, sheetName = "Dendrogram")
insertPlot(DGEwb, sheet = "Dendrogram", width = 2000, height = 3500, units = "px")

lapply(sort(names(finUnList)), function(x){

  # Get DGE results
  res <- finUnList[[x]]$res

  # Add max direction
  res$max_Group <- apply(abs(res[, c("mean_z1", "mean_z2")]), 1, which.max)
  res$dir_Group <- apply(res[, c("mean_z1", "mean_z2")], 1, function(x) sign(x[which.max(abs(x))]))

  # Add reactome pathways
  res$REACTOME_pathways <- gene2Pathway[rownames(res)]

  # Change column names
  colnames(res) <- c("Symbol", "Log2(Fold Change)", "Mean Expression", "Z", "P-Value", "FDR Q-Value", "Beta Value", "Group1 v Vehicle", "Group2 v Vehicle", "Maximum Group", "Group Direction", "REACTOME")
  res$`Group Direction` <- "Up"
  res$`Group Direction`[res$Z < 0] <- "Down"

  # Create simple data frame of chemicals
  chemFram <- data.frame(Groups = c("Group_1", "Group_2"), Chemicals = c(paste(finUnList[[x]]$chems[[1]], collapse = ", "), paste(finUnList[[x]]$chems[[2]], collapse = ", ")))

  # Create worksheet
  addWorksheet(DGEwb, sheetName = paste0(x, "_DGE"))
  writeData(DGEwb, sheet = paste0(x, "_DGE"), chemFram)
  writeData(DGEwb, sheet = paste0(x, "_DGE"), res, startRow = 6)
})

# Add hyperenrichment
lapply(sort(names(finUnList)), function(x){

  # Create simple data frame of chemicals
  chemFram <- data.frame(Groups = c("Group_1", "Group_2"), Chemicals = c(paste(finUnList[[x]]$chems[[1]], collapse = ", "), paste(finUnList[[x]]$chems[[2]], collapse = ", ")))

  HypSub <- HyperEnr[grepl(x, HyperEnr$set),]
  HypSub <- HypSub[HypSub$fdr < 0.5,]
  colnames(HypSub) <- c("set", "P-Value", "FDR Q-Value", "Intersect Size", "DGE Set Size", "Gene Set Size", "Total Genes", "Gene Set Name", "Gene Hits")

  # Make more informative columns
  HypSub$Split <- substr(HypSub$set, 1, 1)
  HypSub$Group <- substr(HypSub$set, nchar(HypSub$set), nchar(HypSub$set))
  HypSub$Direction <- "Up"
  HypSub$Direction[grepl("down", HypSub$set)] <- "Down"
  HypSub <- HypSub[,c("Group", "Direction", "P-Value", "FDR Q-Value", "Gene Set Name", "Intersect Size", "DGE Set Size", "Gene Set Size", "Total Genes", "Gene Hits")]


  addWorksheet(DGEwb, sheetName = paste0(x, "_Pathways"))
  writeData(DGEwb, sheet = paste0(x, "_Pathways"), chemFram)
  writeData(DGEwb, sheet = paste0(x, "_Pathways"), HypSub, startRow = 6)
})
worksheetOrder(DGEwb) <- c(1, unlist(lapply(sort(names(finUnList)), function(x) which(grepl(x, sheets(DGEwb))))))
saveWorkbook(DGEwb, SplitSummaryTableOut, overwrite = T)

# Save RData object of results
save(eSet, finUnList, HyperEnr, gene2Pathway, vNetOut, REAC, dClust, file = rData893a)

# Empty Environment
rm(list = ls())
