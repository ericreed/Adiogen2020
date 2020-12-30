require(edgeR)
require(pheatmap)
require(sva)
require(grid)
require(Biobase)
require(ggplot2)
require(RColorBrewer)
source("parameters.R")

# Set script tag, so we know what script outputs came from
sTag <- "011a"

# Begin

## Read in expression set
eSet <- readRDS(eSetAFile010a)

# Filter for very lowly expressed genes
eSetRem <- eSet[rowMeans(Biobase::exprs(eSet)) >= 1,]

### Create a log2CPM transformed eSet
### edgeR normalize counts
eSetCPM <- eSetRem
dds <- DGEList(Biobase::exprs(eSetCPM), group = eSet$Condition)
dds <- calcNormFactors(dds)
Biobase::exprs(eSetCPM) <- edgeR::cpm(dds, log = TRUE)


### Plot MADvMedianExpression
madVec <- apply(Biobase::exprs(eSetCPM), 1, mad)
medVec <- apply(Biobase::exprs(eSetCPM), 1, median)
madVmed <- data.frame(MADscore = madVec, MedExp = medVec)
png(MADvExpFile011a, 400, 400)
ggplot(madVmed, aes(MedExp, MADscore)) + geom_point()
dev.off()

### Batch correct with combat
eSetCombat <- eSetCPM
Biobase::exprs(eSetCombat) <- ComBat(Biobase::exprs(eSetCombat), batch = eSetCombat$Plate)
clustCom <- createClusterPlots(eSetCombat, 3000, filePre = combatFiles011a, sTag = sTag, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition", "Plate", "Experiment"))

## Plot branches and positions
clustO <- clustCom$tree_col
png(Cluster_Prefilter1, 10000, 800)
plot(clustO)
axis(1, at=1:length(clustCom$tree_col$labels), labels = 1:length(clustCom$tree_col$labels))
dev.off()

## Remove outlier samples (low expression)
pD <- pData(eSet)[colnames(eSet) %in% clustO$labels[clustO$order][190:219], c("Chemical", "Plate")]; write.csv(pD[order(pD$Chemical),], Cluster_RemoveSamps1)
eSet <- eSet[, !colnames(eSet) %in% clustO$labels[clustO$order][190:219] ]

# Filter for very lowly expressed genes
eSetRem <- eSet[rowMeans(Biobase::exprs(eSet)) >= 1,]

### Create a log2CPM transformed eSet
### edgeR normalize counts
eSetCPM <- eSetRem
dds <- DGEList(Biobase::exprs(eSetCPM), group = eSet$Condition)
dds <- calcNormFactors(dds)
Biobase::exprs(eSetCPM) <- cpm(dds, log = TRUE)


### Plot MADvMedianExpression
madVec <- apply(Biobase::exprs(eSetCPM), 1, mad)
medVec <- apply(Biobase::exprs(eSetCPM), 1, median)
madVmed <- data.frame(MADscore = madVec, MedExp = medVec)
png(MADvExpFile_2_011a, 400, 400)
ggplot(madVmed, aes(MedExp, MADscore)) + geom_point()
dev.off()

### Batch correct with combat
eSetCombat <- eSetCPM
Biobase::exprs(eSetCombat) <- ComBat(Biobase::exprs(eSetCombat), batch = eSetCombat$Plate)
clustCom <- createClusterPlots(eSetCombat, 3000, filePre = combatFiles_2_011a, sTag = sTag, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition", "Plate", "Experiment"))

## Plot branches and positions
clustO <- clustCom$tree_col
png(Cluster_Prefilter2, 10000, 800)
plot(clustO)
axis(1, at=1:length(clustCom$tree_col$labels), labels = 1:length(clustCom$tree_col$labels))
dev.off()

## Remove outlier samples (low expression)
pD <- pData(eSet)[colnames(eSet) %in% clustO$labels[clustO$order][190:210], c("Chemical", "Plate")]; write.csv(pD[order(pD$Chemical),], Cluster_RemoveSamps2)
eSet <- eSet[, !colnames(eSet) %in% clustO$labels[clustO$order][190:210] ]

# Filter for very lowly expressed genes
eSetRem <- eSet[rowMeans(Biobase::exprs(eSet)) >= 1,]

### Create a log2CPM transformed eSet
### edgeR normalize counts
eSetCPM <- eSetRem
dds <- DGEList(Biobase::exprs(eSetCPM), group = eSet$Condition)
dds <- calcNormFactors(dds)
Biobase::exprs(eSetCPM) <- cpm(dds, log = TRUE)


### Plot MADvMedianExpression
madVec <- apply(Biobase::exprs(eSetCPM), 1, mad)
medVec <- apply(Biobase::exprs(eSetCPM), 1, median)
madVmed <- data.frame(MADscore = madVec, MedExp = medVec)
png(MADvExpFile_3_011a, 400, 400)
ggplot(madVmed, aes(MedExp, MADscore)) + geom_point()
dev.off()

### Batch correct with combat
eSetCombat <- eSetCPM
Biobase::exprs(eSetCombat) <- ComBat(Biobase::exprs(eSetCombat), batch = eSetCombat$Plate)
clustCom <- createClusterPlots(eSetCombat, 3000, filePre = combatFiles_3_011a, sTag = sTag, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition", "Plate", "Experiment"))

## Plot branches and positions
clustO <- clustCom$tree_col
png(Cluster_Prefilter3, 10000, 800)
plot(clustO)
axis(1, at=1:length(clustCom$tree_col$labels), labels = 1:length(clustCom$tree_col$labels))
dev.off()

## Remove outlier samples (low expression)
pD <- pData(eSet)[colnames(eSet) %in% clustO$labels[clustO$order][17:28], c("Chemical", "Plate")]; write.csv(pD[order(pD$Chemical),], Cluster_RemoveSamps3)
eSet <- eSet[, !colnames(eSet) %in% clustO$labels[clustO$order][17:28] ]

# Filter for very lowly expressed genes
eSetRem <- eSet[rowMeans(Biobase::exprs(eSet)) >= 1,]

### Create a log2CPM transformed eSet
### edgeR normalize counts
eSetCPM <- eSetRem
dds <- DGEList(Biobase::exprs(eSetCPM), group = eSet$Condition)
dds <- calcNormFactors(dds)
Biobase::exprs(eSetCPM) <- cpm(dds, log = TRUE)


### Plot MADvMedianExpression
madVec <- apply(Biobase::exprs(eSetCPM), 1, mad)
medVec <- apply(Biobase::exprs(eSetCPM), 1, median)
madVmed <- data.frame(MADscore = madVec, MedExp = medVec)
png(MADvExpFile_4_011a, 400, 400)
ggplot(madVmed, aes(MedExp, MADscore)) + geom_point()
dev.off()

### Batch correct with combat
eSetCombat <- eSetCPM
Biobase::exprs(eSetCombat) <- ComBat(Biobase::exprs(eSetCombat), batch = eSetCombat$Plate)
clustCom <- createClusterPlots(eSetCombat, 3000, filePre = combatFiles_4_011a, sTag = sTag, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition", "Plate", "Experiment"))

## Plot branches and positions
clustO <- clustCom$tree_col
png(Cluster_Prefilter4, 10000, 800)
plot(clustO)
axis(1, at=1:length(clustCom$tree_col$labels), labels = 1:length(clustCom$tree_col$labels))
dev.off()

# Save
saveRDS(eSetCombat, eSetCombatFile011a)

# Empty Environment
rm(list = ls())
