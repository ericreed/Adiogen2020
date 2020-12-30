require(pROC)
require(openxlsx)
require(ggplot2)
require(Biobase)
source("parameters.R")

# Read in package functions
thisDir <- getwd()
setwd(packDir)
source("R/loadFunctions.R")
setwd(thisDir)

# Set script tag
sTag <- "779a"

# Begin

# Read in file to update chemical names
update <- read.csv(updateFile, stringsAsFactors = F)
update$Abbreviation <- toupper(update$Abbreviation)
update$Abbreviation <- sub("-", "_", update$Abbreviation)

# Read in expression set
eSet <- readRDS(eSetCombatFile011a)

# Remove naive controls
eSet <- eSet[,eSet$Chemical != "NAIVE"]

# Rename vehicle
eSet$PPARg_Mod[grepl("VH", eSet$Chemical)] <- "No"
eSet$Chemical[grepl("VH", eSet$Chemical)] <- "VEHICLE"

# Get features
features <- t(Biobase::exprs(eSet))

# Get pData
pD <- pData(eSet)[, c("Chemical", "PPARg_Mod")]

# Get labelled data
pLab <- pD[pD$PPARg_Mod != "Suspected",]
fLab <- features[rownames(pLab),]
gLab <- pLab$Chemical
cLab <- factor(pLab$PPARg_Mod, levels = c("No", "Yes"))

# Create list with the right objects
labList <- list(
  features = fLab,
  groups = gLab,
  class = cLab
)
names(labList$groups) <- names(labList$class) <- rownames(labList$features)

# Run full random forest

# Filter based on gene expression
runFTest <- function(features, group, intercept = "VEHICLE"){

  # Unique groups
  groups <- unique(group)

  # Set 0 factors to VEHICLE
  group <- factor(group, levels = c(intercept, groups[groups!=intercept]))

  # X matrices for full and reduced model
  Xfull <- model.matrix(~ group)
  Xred <- model.matrix(~ 1, data = as.data.frame(features))

  # Y estimates for full and reduced model
  Bhat <- solve( t(Xfull) %*% Xfull ) %*% t(Xfull) %*% features
  Bbat <- solve( t(Xred) %*% Xred ) %*% t(Xred) %*% features

  Yhat <- Xfull %*% Bhat
  Ybar <- Xred %*% Bbat

  # Y errors for full and reduced model
  SSEfull <- colSums((Yhat-features)^2)
  SSEred <- colSums((Ybar-features)^2)

  # Degrees of freedeom
  nDF <- ncol(Xfull) - ncol(Xred)
  dDF <- (nrow(Xfull)-ncol(Xfull) - 1)

  # Get F-statistics
  Fstat <- ((SSEred - SSEfull)/nDF)/(SSEfull/dDF)

  # Get nominal P and FDR
  DFp <- pf(Fstat, nDF, dDF, lower.tail = F)
  DFfdr <- p.adjust(DFp, method = "BH")

  # Create data.frame of output
  DFout <- data.frame(t(Bhat), Fstat = Fstat, Pval = DFp, FDR = DFfdr)

  return(DFout)
}

Fres <- runFTest(labList$features, labList$group)

# FDR significant genes
FresKeep <- rownames(Fres)[Fres$FDR < 0.05]

# Format data from of filtering results
FresSub <- Fres[FresKeep,]; FresSub <- FresSub[order(FresSub$Pval),]; FresSub <- FresSub[,-1]
colnames(FresSub) <- sub("group", "", colnames(FresSub))
colnames(FresSub)[!colnames(FresSub) %in% c("Fstat", "Pval", "FDR")] <- sapply(colnames(FresSub)[!colnames(FresSub) %in% c("Fstat", "Pval", "FDR")], function(nam){
  update$Update[update$Abbreviation == nam]
})
colnames(FresSub)[!colnames(FresSub) %in% c("Fstat", "Pval", "FDR")]  <- paste0("log2fc", colnames(FresSub)[!colnames(FresSub) %in% c("Fstat", "Pval", "FDR")] )

# Subset lablist
labList$features <- labList$features[,FresKeep]

# Create random forest
RNGkind("L'Ecuyer-CMRG")
set.seed(111)
forestTrain <- build_forest(labList$features,
                            labList$group,
                            labList$class,
                            max.depth = Inf,
                            featPtree = sqrt(ncol(labList$features)),
                            nTrees = 2000,
                            cores = 4,
                            mergeFunction = groupMeans,
                            bagFunction = bag_tree)

# Save
saveRDS(forestTrain, ForestOut779a)

# Get gene importance
gImp <- GiniImportance(forestTrain)
gImpFram <- data.frame(Importance = gImp, row.names = names(gImp))

# Get training threshold
trainEst <- as.data.frame(predict_forest(labList$features,
                                         labList$group, forestTrain, mergeFunction = groupMeans))
# Merge with class info
colnames(trainEst) <- "votes"
trainEst$Chemical <- rownames(trainEst)
info <- unique(data.frame(Chemical = labList$groups, class = labList$class))
trainEst <- merge(trainEst, info)
trainEst <- trainEst[order(trainEst$votes, decreasing = T),]

# Get training threshold
meanYes <- mean(trainEst$votes[trainEst$class == "Yes"]); sdYes <- sd(trainEst$votes[trainEst$class == "Yes"])
meanNo <- mean(trainEst$votes[trainEst$class == "No"]); sdNo <- sd(trainEst$votes[trainEst$class == "No"])
trainThresh <- (sdYes*meanNo + sdNo*meanYes)/(sdYes + sdNo)

# Get OOB Estimates
oobEst <- oob_estimate(forestTrain, features = labList$features, class = labList$group, group = labList$group, mergeFunction = groupMeans); colnames(oobEst) <- c("votes", "Chemical")

# Merge with class info
oobEst <- merge(oobEst, info)
oobEst <- oobEst[order(oobEst$votes),]

# Get OOB threshold
oobEst$F1 <- sapply(1:nrow(oobEst), function(x){

  # Caulcate f1
  sens <- sum(oobEst$class[x:nrow(oobEst)] == "Yes") / sum(oobEst$class == "Yes")
  prec <- sum(oobEst$class[x:nrow(oobEst)] == "Yes") / length(x:nrow(oobEst))
  f1 <- sum(na.omit(c(2 * (sens*prec)/(sens + prec), 0)), na.rm = F)

  return(f1)
})
oobEst$Chemical <- as.character(oobEst$Chemical)
oobEst$Chemical[oobEst$Chemical != "VEHICLE"] <- sapply(oobEst$Chemical[oobEst$Chemical != "VEHICLE"], function(nam) update$Update[update$Abbreviation == nam])
rownames(oobEst) <- oobEst$Chemical; oobEst <- oobEst[,-1]

# Get oobThresh
oobThresh <- mean(c(oobEst$votes[which.max(oobEst$F1)], oobEst$votes[which.max(oobEst$F1)-1]))

# Get Unlabelled data
# Get labelled data
pUnLab <- pD[pD$PPARg_Mod == "Suspected",]
fUnLab <- features[rownames(pUnLab),]
gUnLab <- pUnLab$Chemical

# Create list with the right objects
UnlabList <- list(
  features = fUnLab,
  groups = gUnLab
)
names(labList$groups) <- rownames(labList$features)

# Find estimates of chemicals
UnLabEst <- as.data.frame(predict_forest(UnlabList$features,
                                         UnlabList$group, forestTrain, mergeFunction = groupMeans))
colnames(UnLabEst) <- "votes"
UnLabEst <- UnLabEst[order(UnLabEst$votes, decreasing = T),, drop = F]
UnLabEst$trainPred <- UnLabEst$votes > trainThresh
UnLabEst$oobPred <- UnLabEst$votes > oobThresh
colnames(UnLabEst) <- c("votes", "Pred(Train Threshold)", "Pred(OOB Threshold)")
UnLabEst <- UnLabEst[,c(1, 3)]
rownames(UnLabEst) <- sapply(rownames(UnLabEst), function(nam) update$Update[update$Abbreviation == nam])

# Write results to excel file

## Performance Summary
wb <- createWorkbook()

## Gene Importance
addWorksheet(wb, "Gene_Filtering")
writeData(wb, "Gene_Filtering", FresSub, rowNames = T)

## Gene Importance
addWorksheet(wb, "Gene_Importance")
writeData(wb, "Gene_Importance", gImpFram, rowNames = T)
plotImportance(gImpFram, n = 25)
insertPlot(wb, "Gene_Importance", width = 5, height = 5, xy = NULL, startRow = 1,
           startCol = 4, fileType = "png", units = "in", dpi = 300)
plotImportanceVrank(gImpFram, n = 0) + theme(axis.text = element_text(size = 15))
insertPlot(wb, "Gene_Importance", width = 6, height = 6, xy = NULL, startRow = 50,
           startCol = 4, fileType = "png", units = "in", dpi = 300)

# Thresholds
addWorksheet(wb, "Thresholds_Training")
writeData(wb, "Thresholds_Training", oobEst, rowNames = T)

## Estimates
addWorksheet(wb, "Chemical_Class_Estimates")
writeData(wb,  "Chemical_Class_Estimates", UnLabEst, rowNames = T)
saveWorkbook(wb, excelFileOut779a, overwrite = T)
