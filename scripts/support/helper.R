################### Create Heatmaps
###################################
clusEset <- function(eset, covs = c("PPARg_Mod", "PPARg_Lig", "RXR_Lig", "PPARa_Lig", "Condition"), showgenes = FALSE, cluster_rows=TRUE, cluster_cols = TRUE){
  e <- Biobase::exprs(eset)
  
  #Cluster samples and genes
  if(cluster_cols){
    sDist <- dist(t(e), method = "euclidean")
    sClust <- hclust(sDist, method = "ward.D2")
  } else {sClust = FALSE}
  
  if(cluster_rows){
    gDist <- as.dist(1 - cor(t(e)))
    gClust <- hclust(gDist, method = "ward.D")
  } else {gClust = FALSE}
  
  # Zscore normalize genes
  eS <- t(scale(t(e)))
  eS[eS < -3] <- -3
  eS[eS > 3] <- 3
  
  # Create annotation
  annot <- pData(eset)
  annot <- annot[colnames(e),]
  
  platCols <- list(Plate = brewer.pal(6, "Pastel1"),
                   PPARg_Mod = c(brewer.pal(3, "Set1")[1:3]),
                   PPARg_Lig = "black",
                   RXR_Lig = "black",
                   PPARa_Lig = "black",
                   Condition = c("lightgrey", "lightgreen", brewer.pal(4, "Set1")[4]))
  names(platCols$PPARg_Mod) <- c("Yes", "No", "Suspected")
  names(platCols$PPARg_Lig) <- c("Yes")
  names(platCols$PPARa_Lig) <- c("Yes")
  names(platCols$RXR_Lig) <- c("Yes")
  names(platCols$Plate) <- paste0("plate", 1:6)
  names(platCols$Condition) <- c("Naive Control", "Vehicle Control", "Treated")
  
  #can do show_rownames=FALSE OR TRUE
  pheatmap(eS, color = rev(brewer.pal(11, "RdBu")), cluster_rows = gClust, cluster_cols = sClust, annotation_col = annot[,covs], show_rownames = showgenes, labels_col = annot$Chemical, annotation_colors = platCols)
}

################### MAD Score
##############################
madFilter <- function(eset, filter = 3000){
  madVec <- apply(Biobase::exprs(eset), 1, mad)
  madVec <- sort(madVec, decreasing = TRUE)
  eSetMAD <- eset[names(madVec)[1:filter],]
  return(eSetMAD)
}

################### Expression
##############################
expFilter <- function(eset, filter = 3000, func = median){
  medVec <- apply(Biobase::exprs(eset), 1, func)
  medVec <- sort(medVec, decreasing = TRUE)
  eSetMED <- eset[names(medVec)[1:filter],]
  return(eSetMED)
}

################### Cluster Plots
#################################

createClusterPlots <- function(eset, nGenes, covs = c("Condition", "Obesogen", "PPARg", "Plate"), filePre, sTag){
  # Create file names to write out
  clustMADplotFile <- paste0(filePre, "clustMADheatmap_", nGenes, "genes_", sTag, ".png")
  clustExpplotFile <- paste0(filePre, "clustExpheatmap_", nGenes, "genes_", sTag, ".png")
  pcaMADplotFile <- paste0(filePre, "pcaMADplot_", nGenes, "genes_", sTag, ".png")
  pcaExpplotFile <- paste0(filePre, "pcaExpplot_", nGenes, "genes_", sTag, ".png")
  
  ### Filter for top 3K MAD scores and cluster
  eSetMAD <- madFilter(eset, filter = nGenes)
  #### Create PCA plot
  png(pcaMADplotFile, width = 500, height = 500)
  grid.draw(pcaPlot(eSetMAD))
  dev.off()
  ### Cluster top 3K mad genes
  png(clustMADplotFile, width = 3000, height = 1000)
  clusEset(eSetMAD, covs = covs)
  dev.off()
  
  ### Filter for top 3K MAD Expressed genes and cluster
  eSetMED <- expFilter(eset, filter = nGenes)
  ### Create PCA plot
  png(pcaExpplotFile, width = 500, height = 500)
  grid.draw(pcaPlot(eSetMED))
  dev.off()
  ### Cluster top 3K mad genes
  png(clustExpplotFile, width = 3000, height = 1000)
  out <- clusEset(eSetMED, covs = covs)
  dev.off()
  return(out)
}

# Run LIMMA
###########

runLIMMA <- function(eSet, control, outDir, sTag, outGSfile, cores = 1, ran = "Plate"){
  chemT <- unique(eSet$Chemical[eSet$Condition!=control])
  cOut <- mclapply(chemT, function(x){
    
    # Run Limma
    eSetSub <- eSet[,eSet$Chemical == x | eSet$Condition == control]
    e <- Biobase::exprs(eSetSub)
    eSetSub$Condition[eSetSub$Condition!=control] <- "Treated"
    eSetSub$Condition <- factor(eSetSub$Condition, levels = c(control, "Treated"))
    design <- model.matrix(~ Condition, pData(eSetSub))
    colnames(design)[ncol(design)] <- "Treated"
    # Account for Tech Reps
    corfit <- duplicateCorrelation(e,design,block=pData(eSetSub)[,ran])
    fit <- lmFit(e, design, block = pData(eSetSub)[,ran], correlation = corfit$consensus)
    fit <- eBayes(fit, trend=TRUE)
    out <- topTable(fit, coef="Treated", number = Inf, sort.by = "P")
    out$high.class <- c("Vehicle Control", "Treated")[(out$logFC>0)+1]
    out$gene <- rownames(out)
    # Create gene signature
    Signature <- out[, c("gene", "P.Value", "logFC", "high.class", "AveExpr")]
    rownames(Signature) <- NULL
    colnames(Signature) <- c("ID", "P.Value", "Score", "high.class", "baseMean")
    Signature$high.class <- unlist(Signature$high.class)
    MetaData <- data.frame(info =  c("Vehicle Control", "Treated"), stringsAsFactors = FALSE)
    row.names(MetaData) <- c("Control", "Treatment")
    GS <- GeneSig(SigType = "gene_rank", MetaData = MetaData, Signature = Signature)
    
    # Write Gene Signature
    write.GeneSig(GS, paste(outGSfile, x, "Vehicle", sTag, sep = "_"), outDir, writeRDS = TRUE)
    
    return(out)
  }, mc.cores = cores)
  names(cOut) <- chemT
  
  # return names list
  return(cOut)
}

# Write list of results to xlsx files
#####################################
writeXLSX <- function(namedList, fileName){
  if(file.exists(fileName)) file.remove(fileName)
  
  # Write data
  wb <- openxlsx::createWorkbook()
  lapply(1:length(namedList), function(x){
    out <- namedList[[x]]
    nam <- names(namedList)[x]
    
    addWorksheet(wb, nam)
    writeData(wb, sheet = nam, x = out, rowNames = T)
  })
  openxlsx::saveWorkbook(wb, fileName, overwrite = TRUE)
}

# Read in GMT files
getgeneLists <-function(GMTfile, rem = c(1,2)){
  GMTVec <- readLines(GMTfile)
  GMTList <- strsplit(GMTVec, "\t")
  GMTNames <- unlist(lapply(GMTList, function(x) x[1]))
  GeneSets <- lapply(GMTList, function(x) x[-rem])
  names(GeneSets) <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", GMTNames, perl=TRUE)
  return(GeneSets)
}

# Plot library sizes by plate

boxLS <- function(eSet, sepCol = "Plate", ylimits = c(0, 6.5), addPoints = FALSE, Ny = 6.5, sortMean = F, eSetRef = NULL, RefName = NULL, untrans = FALSE, ylab = "Log10( Library Size )"){
  
  # Data
  df <- data.frame(sepCol = pData(eSet)[,sepCol], LS = log(colSums(Biobase::exprs(eSet)), 10))
  df$sepCol <- factor(df$sepCol, levels = unique(df$sepCol))
  
  # Add Refereace
  if(!is.null(eSetRef)){
    dfAdd <- data.frame(sepCol = RefName, LS = log(colSums(Biobase::exprs(eSetRef)), 10))
    df$sepCol <- as.character(df$sepCol)
    df <- rbind(df, dfAdd)
    df$sepCol <- factor(df$sepCol, levels = unique(df$sepCol))
  }
  
  # Sort Mean?
  if(sortMean){
    dfMean <- df %>% group_by(sepCol) %>% summarise(meanLS = mean(LS))
    dfMean <- dfMean[order(dfMean$meanLS),]
    df$sepCol <- factor(df$sepCol, levels = as.character(dfMean$sepCol))
  }
  colnames(df)[1] <- sepCol
  
  # Untransform plot
  if(untrans) df$LS <- 10^df$LS
  
  # Counts
  tab <- table(df[,sepCol])
  tab <- tab[levels(df[,sepCol])]
  dfCounts <-data.frame(N = tab)
  colnames(dfCounts)[1] <- "SepCol"
  dfCounts[,sepCol] <- factor(dfCounts$SepCol, levels = levels(df[,sepCol]))
  
  p <- ggplot(df, aes_string(x = sepCol, y = "LS")) +
    geom_boxplot() +
    geom_text(data = dfCounts, aes(x = 1:length(unique(df[,sepCol])), y = Ny, label = N.Freq)) +
    scale_y_continuous(name = ylab, limits = ylimits) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
      )
  if(addPoints) p <- p + geom_point(alpha = 0.6)
  return(p)
}

# Plot importance vs ranking

plotImportanceVrank <- function(gImpFram, n = 5){
  colnames(gImpFram) <- "Importance"
  gImpFram$Rank <- 1:nrow(gImpFram)
  gImpFram$Gene[1:n] <- rownames(gImpFram)[1:n]
  nameFram <- rownames(gImpFram)
  
  p <- ggplot(gImpFram, aes(x = Rank, y = Importance)) +
    geom_point() + 
    theme_bw()
  return(p)
}

# Wrapper for creating expression sets
to.eSet <- function (mat, pdat, fdat) {
  mat <- as.matrix(mat)
  if (!is.data.frame(pdat)) 
    stop("pdat must be a data frame")
  if (!is.data.frame(fdat)) 
    stop("fdat must be a data frame")
  if (nrow(fdat) != nrow(mat)) 
    stop("nrow(fdat) must equal nrow(mat)")
  if (nrow(pdat) != ncol(mat)) 
    stop("nrow(pdat) must equal ncol(mat)")
  if (!all(rownames(fdat) == rownames(mat))) {
    warning("fdat rownames and mat rownames do not match, setting fdat rownames to mat rownames")
    rownames(fdat) <- rownames(mat)
  }
  if (!all(rownames(pdat) == colnames(mat))) {
    warning("pdat rownames and mat colnames do not match, setting fdat rownames to mat rownames")
    rownames(pdat) <- colnames(mat)
  }
  fMetaData <- data.frame(labelDescription = colnames(fdat), 
                          row.names = colnames(fdat))
  featureData <- new("AnnotatedDataFrame", data = fdat, varMetadata = fMetaData)
  pMetaData <- data.frame(labelDescription = colnames(pdat), 
                          row.names = colnames(pdat))
  phenoData <- new("AnnotatedDataFrame", data = pdat, varMetadata = pMetaData)
  eSet <- ExpressionSet(assayData = mat, featureData = featureData, 
                        phenoData = phenoData, annotation = "")
  return(eSet)
}

################### PCA plot
##############################

pcaPlot <- function(eset) {
  e <- Biobase::exprs(eset)
  p <- princomp(e)
  df <- as.data.frame(p$loadings[,1:2])
  perVar <- round((p$sdev^2/sum(p$sdev^2)*100)[1:2], 2)
  df$Plate <- pData(eset)$Plate
  out <- ggplot(df, aes(Comp.1, Comp.2)) +
    geom_point(aes(color = Plate)) +
    scale_x_continuous(name = paste0("PC1: %Var = ", perVar[1]))+
    scale_y_continuous(name = paste0("PC2: %Var = ", perVar[2]))
  return(out)
}

################### Ordering HClust branches
############################################

hcopt <- function (d, HC = NULL, method = "ward.D", members = NULL) {
  if (is.null(HC)) {
    HC <- hclust(d, method = method, members = members)
  }
  if (length(HC$order) > 2) {
    ORD <- order.optimal(d, merge = HC$merge)
    HC$merge <- ORD$merge
    HC$order <- ORD$order
  }
  HC
}