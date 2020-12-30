require(edgeR)
require(pheatmap)
require(dplyr)
source("parameters.R")

# Set script tag, so we know what script outputs came from
sTag <- "010a"

# Main directories
baseDir <- ".."
datDir <- file.path(baseDir, "data")
resDir <- file.path(baseDir, "results")
plotDir <- file.path(resDir, "plots")
tabDir <- file.path(resDir, "tables")

## Subdirectories
esetDir <- file.path(datDir, "eSets")

# Files to read in
eSetFiles <- list.files(esetDir, pattern = "_3DGE_LERem_plate")

# Files to write out

# Read in eSets
eSetList <- lapply(eSetFiles, function(x) readRDS(file.path(esetDir, x)))
names(eSetList) <- gsub("Obesogen_3DGE_LERem_|_00[[:digit:]].rds", "", eSetFiles)

# Add plateID to each eSet
eSetList <- lapply(1:length(eSetList), function(x){
  pData(eSetList[[x]])$Plate = names(eSetList)[x]
  colnames(eSetList[[x]]) <- paste(colnames(eSetList[[x]]), names(eSetList)[x], sep = "_")
  return(eSetList[[x]])})

# Combine into a single eset
eFram <- do.call(cbind, lapply(eSetList, function(x) exprs(x)))
pFram <- do.call(rbind, lapply(eSetList, function(x) pData(x)))
fFram <- fData(eSetList[[1]])
eSet <- to.eSet(eFram, pFram, fFram)

# Add Plate IDs to vehicles
eSet$Chemical[eSet$Chemical == "VH"] <- paste0(eSet$Chemical[eSet$Chemical == "VH"], 
                                               eSet$Plate[eSet$Chemical == "VH"])

## Write to file
saveRDS(eSet, eSetRawFile010a)

## Cluster the data based on naive normalization (log2CPM)
### Create a log2CPM transformed eSet
eSetCPM <- eSet
exprs(eSetCPM) <- cpm(eSetCPM, log = TRUE)


#### Filer for top 3000 MAD Score genes
eSetCPM_MAD <- madFilter(eSetCPM)
### Create Heatmap
png(HM_CPM_MADFile010a, width = 4000, height = 1000)
outClustMed <- clusEset(eSetCPM_MAD, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition", "Plate", "Experiment"))
dev.off()

#### Filer for top 3000 MEDIAN EXPRESSION genes
eSetCPM_MED <- expFilter(eSetCPM)
### Create Heatmap
png(HM_CPM_MEDFile010a, width = 4000, height = 1000)
outClustMed <- clusEset(eSetCPM_MED, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition", "Plate", "Experiment"))
dev.off()

### Add library size to each sample
eSum <- log(colSums(exprs(eSet)) + 1, 10)
pData(eSet)$Log10LibrarySize <- eSum


# Separate into exp A and exp B esets and write to file
## Exp. A
eSetA <- eSet[,eSet$Experiment == "A"]
saveRDS(eSetA, eSetAFile010a)

## Exp. B
eSetB <- eSet[,eSet$Experiment == "B"]
saveRDS(eSetB, eSetBFile010a)

# Explore Plates Separately

# Plot library sizes

# All chems
## By Chemical
eSetSub <- eSetA
eSetSub$Chemical[grepl("VH", eSetSub$Chemical)] <- "VEHICLE"
tabChem <- table(eSetSub$Chemical); tabChem <- tabChem[tabChem>2]; eSetSub <- eSetSub[,eSetSub$Chemical %in% names(tabChem)]
LSSub <- boxLS(eSetSub, ylimits = c(1, 7), addPoints = T, Ny = 7, sepCol = "Chemical", sortMean = T)
png(LibrarySizeByChem, width = 1200, height = 500)
LSSub
dev.off()

# Vehicles Only
## By Plate
eSetAv <- eSetA[, grepl("VH", eSetA$Chemical)]
LSAv <- boxLS(eSetAv, ylimits = c(3, 6.2), addPoints = T, Ny = 6)
png(LibrarySizeVehicles, width = 500, height = 500)
LSAv
dev.off()

# Plate 5
eSetP5 <- eSetA[, eSetA$Plate == "plate5" ]
eSetP5$Chemical[grepl("VH", eSetP5$Chemical)] <- "VEHICLE"
tabChem <- table(eSetP5$Chemical); tabChem <- tabChem[tabChem>2]; eSetP5 <- eSetP5[,eSetP5$Chemical %in% names(tabChem)]
LSP5 <- boxLS(eSetP5, ylimits = c(2, 7), addPoints = T, Ny = 7, sepCol = "Chemical", sortMean = T)
png(LibrarySizePlate5, width = 900, height = 500)
LSP5
dev.off()

# Plate 6
eSetP6 <- eSetA[, eSetA$Plate == "plate6" ]
eSetP6$Chemical[grepl("VH", eSetP6$Chemical)] <- "VEHICLE"
tabChem <- table(eSetP6$Chemical); tabChem <- tabChem[tabChem>2]; eSetP6 <- eSetP6[,eSetP6$Chemical %in% names(tabChem)]
LSP6 <- boxLS(eSetP6, ylimits = c(1, 7), addPoints = T, Ny = 7, sepCol = "Chemical", sortMean = T)
png(LibrarySizePlate6, width = 500, height = 500)
LSP6
dev.off()

# Empty Environment
rm(list = ls())
