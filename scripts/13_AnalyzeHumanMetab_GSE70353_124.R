source("parameters.R")

require(Biobase)
require(GSVA)
require(ggplot2)
require(reshape2)
require(limma)
require(openxlsx)
require(pheatmap)
require(RColorBrewer)
require(ggdendro)
require(dendextend)
require(ppcor)
require(dplyr)
require(grid)
require(gridExtra)
require(png)

# Read in tax results
load(rData893a)

# Read in variable names
variableNames <- read.csv(varInfoFile, row.names = 1)
rownames(variableNames) <- toupper(rownames(variableNames))

# Get gene lists
gLists <- unlist(lapply(finUnList, function(x){
  x <- x$res
  x <- x[x$adj.P.Val < 0.1 & abs(x$logFC) > 1,]
  x$diff <- apply(abs(x[, c("mean_z1", "mean_z2")]), 1, which.max)
  list(up1 = rownames(x)[x$diff == 1 & x$mean_z1 > 0],
       down1 = rownames(x)[x$diff == 1 & x$mean_z1 < 0],
       up2 = rownames(x)[x$diff == 2 & x$mean_z2 > 0],
       down2 = rownames(x)[x$diff == 2 & x$mean_z2 < 0])
}), recursive = F)
names(gLists) <- sub("[.]", "_", names(gLists))
names(gLists) <- sub("up", "up_group", names(gLists))
names(gLists) <- sub("down", "down_group", names(gLists))

# Read in expression values
eSet <- readRDS(eSet_Human_124)
rownames(eSet) <- toupper(rownames(eSet))

# Perform ssGSEA
if(!file.exists(eSet_Human_ssGSEA_tax_124)) {
  gseSet <- gsva(eSet, gLists, method="gsva")
  saveRDS(gseSet, eSet_Human_ssGSEA_tax_124)
} else {
  gseSet <- readRDS(eSet_Human_ssGSEA_tax_124)
}

# Plot histograms of measurements
pDmelt <- melt(pData(gseSet)[,4:27])

ggplot(pDmelt,aes(x=value)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free")

varUse <- colnames(pData(eSet))[-c(1:3, 5:8, 10, 12, 15, 18, 21, 26, 27)]

# Add partial correlation
design <- pData(gseSet)[,varUse]
res <- do.call(rbind, mclapply(rownames(gseSet), function(gs){

    # Create data.frame
    df <- data.frame(gs = exprs(gseSet)[gs,], design)
    df <- df[complete.cases(df),]
    pcor <- pcor(df, method = "spearman")

    # Perform partial correlation test
    pcorEst <- pcor$estimate["gs",-c(1:2)]
    pcorP <- pcor$p.value["gs",-c(1:2)]

    return(list(estimate = pcorEst,
                p.value = pcorP
    ))
}))

# Get stat matrices and calc FDR
resEstimate <- do.call(rbind, res[,1])
resPvalue <- do.call(rbind, res[,2])
resFDR <- matrix(p.adjust(resPvalue, method = "BH"), ncol = ncol(resPvalue))
rownames(resEstimate) <- rownames(resPvalue) <- rownames(resFDR) <- rownames(gseSet)
colnames(resFDR) <- colnames(resEstimate)

# Create matrix of results
resCombined <- do.call(cbind, lapply(colnames(resEstimate), function(colN) {
    dfSub <- data.frame(
        pcor = resEstimate[,colN],
        pval = resPvalue[,colN],
        fdr = resFDR[,colN]
    )
}))
rowInfo <- strsplit(rownames(resCombined), "_")

# Combine names
resCombined$group <- unlist(lapply(rowInfo, function(x) paste0(x[1], gsub("group", "", x[3]))))
resCombined$direction <- unlist(lapply(rowInfo, function(x) x[2]))
resCombined$genes <- unlist(lapply(rownames(resCombined), function(x) paste(gLists[[x]], collapse = "; ")))

# Write an xlsx file
write.xlsx(resCombined, Human_ssGSEA_tax_res_124, rowNames = T)


scatterGS <- function(gs, marker){

    # Get measured values
    df <- data.frame(marker = pData(gseSet)[,marker], gs = exprs(gseSet)[gs,])

    # Get long title for variable
    namLong <- as.character(variableNames[toupper(marker),1])
    print(namLong)

    png(paste0(Human_ssGSEA_tax_res_prefix, "scatterplot", "_", marker, "_", gs, ".png"), width = 700, height = 700)
    print(ggplot(df, aes(y = gs, x = marker)) +
        theme_bw() +
            geom_point() +
            scale_x_continuous(name = namLong) +
            scale_y_continuous(name = "Enirchment Score") +
            ggtitle(gs) +
            geom_smooth(method = "lm", se = F) +
            theme(
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 20),
                plot.title = element_text(size = 23)
            )
        )
    dev.off()
}

scatterGS("A_up_group2", "p_adipon")
scatterGS("A_down_group2", "p_adipon")
scatterGS("J_up_group2", "p_adipon")
scatterGS("J_down_group2", "p_adipon")
scatterGS("H_up_group2", "p_adipon")

# Get coordinates of the dendrogram
dendroCoords <- unique(dendro_data(dClust)$segments)
dendroCoords <- dendroCoords[order(dendroCoords$x),]
dendroCoords <- dendroCoords[order(dendroCoords$y, decreasing = T),]
dendroCoords <- dendroCoords[dendroCoords$yend != 0 & dendroCoords$y != 1,]

# Get vertical and horizontal lines
dendroVert <- unique(dendroCoords[dendroCoords$x == dendroCoords$xend, c("x", "y", "yend")])

# Add split and group IDs
dendroVert$split <- rep(names(finUnList), each = 2)
dendroVert$group <- rep(c(1, 2) , length(names(finUnList)))

# Add cordinates for adding rectangles
dendroVert$x1 <- dendroVert$x - 0.4
dendroVert$x2 <- dendroVert$x + 0.4
dendroVert$y1 <- dendroVert$y - 0.2
dendroVert$y2 <- dendroVert$y + 0.2

# Reformat to draw polygons (will be flipped along X axis)
triUp <- do.call(rbind, lapply(1:nrow(dendroVert), function(ind) {

  # Get row
  dRow <- dendroVert[ind,]
  dRow$direction <- "up"
  dRow$gTri <- paste(dRow$split, dRow$group, dRow$direction, sep = "_")

  # Replicate 3 times
  dRow4 <- do.call(rbind, lapply(1:3, function(x) dRow))

  # Add coordinates
  dRow4$xT1 <- c(dRow$x2, dRow$x2, dRow$x1)
  dRow4$yT1 <- c(dRow$y1, dRow$y2, dRow$y2)

  return(dRow4)

}))

# Reformat to draw polygons
triDown <- do.call(rbind, lapply(1:nrow(dendroVert), function(ind) {

  # Get row
  dRow <- dendroVert[ind,]
  dRow$direction <- "down"
  dRow$gTri <- paste(dRow$split, dRow$group, dRow$direction, sep = "_")

  # Replicate 3 times
  dRow4 <- do.call(rbind, lapply(1:3, function(x) dRow))

  # Add coordinates
  dRow4$xT1 <- c(dRow$x2, dRow$x1, dRow$x1)
  dRow4$yT1 <- c(dRow$y1, dRow$y1, dRow$y2)

  return(dRow4)

}))

# Combine
dendTri <- rbind(triUp, triDown)

# Add Labels of leaves
leaf_labels = get_leaves_attr(dClust, "label")

# Calculate FDR thresholds
#FDR25 <- min(abs(resEstimate)[resFDR < 0.25]) - 1E-6
FDR10 <- min(abs(resEstimate)[resFDR < 0.1]) - 1E-6
FDR05 <- min(abs(resEstimate)[resFDR < 0.05]) - 1E-6
FDR01 <- min(abs(resEstimate)[resFDR < 0.01]) - 1E-6

# Create plot

# Get rownames to plot
resNames <- colnames(resEstimate)

# Read in key plot
keyIMG <- readPNG(Human_key_file)

# Plot the dendrogram
resList <- lapply(resNames, function(resName){

  # Get Results
  res <- as.data.frame(resEstimate[,resName,drop = F])
  colnames(res) <- "pcor"

  #res$pCutoff <- - FDR05
  #changeMat <- rbind(c(-FDR05, -FDR10, -FDR10), c(-FDR10, -FDR25, -FDR25), c(-FDR25, FDR25, 0), c(FDR25, FDR10, FDR25), c(FDR10, FDR05, FDR10), c(FDR05, Inf, FDR05))

  res$pCutoff <- - FDR01
  changeMat <- rbind(c(-FDR01, -FDR05, -FDR05), c(-FDR05, -FDR10, -FDR10), c(-FDR10, FDR10, 0), c(FDR10, FDR05, FDR10), c(FDR05, FDR01, FDR05), c(FDR01, Inf, FDR01))


  for(i in seq(nrow(changeMat))) res$pCutoff[res$pcor > changeMat[i,1] & res$pcor < changeMat[i,2]] <- changeMat[i,3]

  # Get split information
  rowSplit <- strsplit(rownames(res), "_")
  res$split <- unlist(lapply(rowSplit, function(x) x[1]))
  res$direction <- unlist(lapply(rowSplit, function(x) x[2]))
  res$group <- unlist(lapply(rowSplit, function(x) gsub("group", "", x[3])))

  # Create splitGroup variable
  res$splitgroup <- paste0(res$split, res$group)

  # Merge plots with the triangle
  # resTri <- merge(dendTri, res, all.x = TRUE)
  resTri <- merge(dendTri, res)

  dend_segs <- dendro_data(dClust, type = "rectangle")$segments
  dend_horiz <- unique(dend_segs[dend_segs$y == dend_segs$yend & dend_segs$y > 1, c("x", "y")])
  dend_horiz <- dend_horiz[order(dend_horiz$x),]
  dend_horiz <- dend_horiz[order(dend_horiz$y, decreasing = T),]
  dend_horiz$label <- names(finUnList)

  # Get title
  varNameUse <- as.character(variableNames[toupper(resName),1])
  TITLE <- paste0(varNameUse, " gene set projection")

  # Get number of leaves
  nLeaves <- attributes(dClust)$members

  # Create dendrogram
  pdendro <- ggdendrogram(dClust) +
    geom_segment(data = resTri, aes(x = x, xend = x, y = y, yend = 0), linetype = "dotted") +
    geom_label(data = dend_horiz, aes(x = x, y = y, label = label), size = 6, color = "red") +
    geom_polygon(data = resTri, aes(x = xT1, y = yT1, group = gTri, fill = pCutoff), color = "black") +
    geom_text(data = resTri, aes(x = x2, y = y2, label = splitgroup), size = 6, nudge_y = 0.11, nudge_x = 0.13) +
    scale_x_reverse(breaks = 1:length(leaf_labels), label = leaf_labels, position = "bottom") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name= "Correlation") +
    geom_polygon(x = c(42.4, 42.4, 41.6), y = c(8.8, 9.2, 8.8), group = c(1, 1, 1), fill = "red", color = "black") +
    ggtitle(TITLE) +
    theme(
      legend.position="none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(vjust = 0, size = 20),
      legend.title.align = 0.5,
      plot.title = element_text(size=22),
      legend.text=element_text(size=13)
    )

    # Plot the values as a ranking
    # Subset for unique res
    resUnique <- unique(resTri[,c("splitgroup", "direction", "x", "pcor", "pCutoff")])

    # Up on top
    resUnique$direction <- toupper(resUnique$direction)
    resUnique$direction <- factor(resUnique$direction, levels = c("UP", "DOWN"))

    # Jitter x-coord
    resUniqueNoDirection <- unique(resUnique[,c("splitgroup", "x")])
    resUniqueNoDirection$xJitter <- resUniqueNoDirection$x
    resUniqueNoDirection$xJitter[duplicated(resUniqueNoDirection$xJitter)] <- resUniqueNoDirection$xJitter[duplicated(resUniqueNoDirection$xJitter)] + 0.5
    resUnique <- merge(resUnique, resUniqueNoDirection)

    # Set order to x coord
    resUnique <- resUnique[order(resUnique$xJitter),]
    resUnique$splitgroup <- factor(resUnique$splitgroup, levels = unique(resUnique$splitgroup))

    # Spread 'em
    pscatter <- ggplot(resUnique, aes(x = xJitter, y = pcor)) +
      geom_hline(yintercept=0, alpha = 0.9) +
      geom_hline(yintercept=FDR10, alpha = 0.5) +
      geom_hline(yintercept=-FDR10, alpha = 0.5) +
      geom_point(aes(colour= pCutoff), size = 8) +
      geom_point(shape = 1, size = 6, colour = "black") +
      scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      facet_grid(direction~.) +
      scale_y_continuous(name = "Partial Correlation", breaks = c(-.2, -.15, -.1, 0, .1, .15, .2)) +
      scale_x_reverse(expand = c(0,0), breaks = unique(resUnique$xJitter), labels = levels(resUnique$splitgroup), limits = c(nLeaves, 0)) +
      theme_bw() +
      theme(
        panel.spacing = unit(2, "lines"),
        legend.position = "none",
        panel.grid.major.x = element_line(colour = "black"),
        strip.text.y = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.title.align = 0.5,
        panel.border = element_rect(colour = "black", fill=NA, size=2)
      )

  fileSuffix <- resName

  # Write to file
  plotFile <- paste0(Human_ssGSEA_tax_res_prefix, fileSuffix, ".png")
  png(plotFile , width = 1200, height = 1000)
  grid.arrange(pdendro, pscatter, layout_matrix = matrix(c(1, 1, 2), ncol = 1))
  dev.off()

  # Read in plot
  plotIMG <- readPNG(plotFile)
  resDim <- dim(plotIMG)[2:1]

  png(plotFile , width = 1200, height = 1000)
  plot(1,1,xlim=c(1,resDim[1]),ylim=c(1,resDim[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(plotIMG,1,1,resDim[1],resDim[2])
  rasterImage(keyIMG,1032,725,1241,944)
  dev.off()

  # Subset res for output
  resUnique$varName <- varNameUse

  return(resUnique)

})
names(resList) <- resNames

# Create a cool plot of all of the Results

# # Get column order (Adiponectin order)
colOrder <- rev(levels(resList$p_adipon$splitgroup))

# Get a frame of all results
resFrame <- do.call(rbind, resList)

# Remove instances where there are no significant histograms
resKEEP <- table(resFrame$varName, abs(resFrame$pCutoff) >= FDR10 )[,"TRUE"] > 0
resKEEP <- names(resKEEP)[resKEEP]

resFrame <- resFrame[resFrame$varName %in% resKEEP,]

# Get just J2
J2 <- resFrame[resFrame$splitgroup == "J2",]

# Get row order based on J2 results and for Adipnectin to top
J2$pcor[J2$direction == "DOWN"] <- J2$pcor[J2$direction == "DOWN"] * -1
J2Means <- J2 %>% group_by(varName) %>% summarise(pcor = sum(pcor))
J2Means <- J2Means[order(J2Means$pcor, decreasing = T),]
J2Means <- J2Means[c(2, 1, 3:nrow(J2Means)),]

# Set factors
resFrame$splitgroup <- factor(resFrame$splitgroup, levels = colOrder)
resFrame$varName <- factor(resFrame$varName, levels = rev(J2Means$varName))

# Create triangles
resFrame$x <- as.numeric(resFrame$splitgroup)
resFrame$y <- as.numeric(resFrame$varName)

# Triple it
resFramTri <- do.call(rbind, lapply(seq(nrow(resFrame)), function(nr){

    # Get row
    resRow <- resFrame[nr,]
    resRow$splitgroupdirection <- paste0(resRow$splitgroup, resRow$direction, resRow$varName)

    # Get direction
    up <- c(-0.45, 0.45)[as.numeric(resRow$direction == "UP") + 1]

    # Set trianle coordinates
    dfAdd <- data.frame(xT1 = c(resRow$x - up, resRow$x - up, resRow$x + up),
                        yT1 = c(resRow$y + up, resRow$y - up, resRow$y + up),
                        splitgroupdirection = resRow$splitgroupdirection)

    # Merge
    resRows <- merge(dfAdd, resRow)

    return(resRows)

}))

resFramTri <- resFramTri

# draw plot
png(Human_ssGSEA_Heatmap, width = 1300, height = 310)
print(ggplot(data = resFramTri) +
    geom_polygon(data = resFramTri,
        aes(x = xT1, y = yT1, group = splitgroupdirection, fill = pCutoff),
        color = "black", size = 0.2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name= "Correlation") +
    scale_x_continuous(limits = c(0.5,27.5), expand = c(0, 0), breaks = 1:length(levels(resFrame$splitgroup)), labels = levels(resFrame$splitgroup), position = "top") +
    scale_y_continuous(limits = c(0.5, length(resKEEP) + 0.5), expand = c(0, 0), breaks = 1:length(levels(resFrame$varName)), labels = levels(resFrame$varName), position = "right") +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid.major = element_line(colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(vjust = 0, size = 18),
        axis.text.y = element_text(size = 15),
        axis.ticks.y = element_blank()
    ))
dev.off()

# Create a connector segment plot
resSeg <- unique(resList[[1]][,c("splitgroup", "xJitter", "x")])
resSeg$x2 <- as.numeric(resSeg$splitgroup)*1.5+6
resSeg$y1 <- 1
resSeg$y2 <- 0
colnames(resSeg) <- c("splitgroup", "x1", "x", "x2", "y1", "y2")

png(Human_ssGSEA_Lines, height = 100, width = 1300)
ggplot(data = resSeg) +
    scale_x_reverse() +
    geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2)) +
    theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
dev.off()

png(Human_ssGSEA_Lines2, height = 50, width = 1300)
ggplot(data = resSeg) +
    scale_x_reverse() +
    geom_segment(aes(x = x, xend = x1, y = y1, yend = y2), linetype = "dotted", size = 0.87) +
    theme(
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
dev.off()
