# Create dendrogram
mat2dendro <- function(finUnList, splitNames = LETTERS[1:length(finUnList)]){


  # Create Matric of results
  mat <- matrix(0, nrow = length(finUnList), ncol = ncol(Dfram))
  colnames(mat) <- colnames(Dfram)
  for(i in 1:length(finUnList)){
    mat[i,finUnList[[i]]$chems[[1]]] <- 1
    mat[i,finUnList[[i]]$chems[[2]]] <- 2
  }

  ## collapse mat
  matCollapse <- sort(apply(mat, 2, function(x) paste(x[x!=0], collapse = "")))
  matUnique <- unique(matCollapse)

  # Get branchlist
  bList <- c()
  j <- 1
  for(i in 1:max(nchar(matUnique))){
    sLength <- matUnique[nchar(matUnique) >= i]
    sLength <- unique(substr(sLength, 1, i))
    for(k in sLength){bList[j] <- k; j <- j+1}
  }

  # Add edgenames
  for(i in 2:nrow(mat)){
    matNow <-  mat[i,]
    wn0 <- which(matNow!=0)[1]
    matPrev <- mat[1:(i-1),wn0]
    matPrev <- paste(matPrev[matPrev!=0], collapse = "")
    names(splitNames)[i] <- matPrev
  }


  # Add list
  val <- max(nchar(matUnique)) + 1
  aList <- list()
  attributes(aList)<-list(members=length(matCollapse),height=val, midpoint = (length(matCollapse)-1)/2, label =splitNames[1])

  for(i in bList){
    # Add element
    iSplit <- unlist(strsplit(i, ""))
    iPaste <- paste0("aList", paste(paste0("[[", iSplit, "]]"), collapse = ""))
    eval(parse(text=
                 paste0(iPaste, "<-list()")
    ))

    # Add attributes
    members <- sum(substr(matCollapse, 1, nchar(i)) == i)
    height <- val-nchar(i)
    if(!any(substr(bList, 1, nchar(bList)-1) == i)) height = 1
    att <- list(members = members, height = height, midpoint = (members-1)/2, label = splitNames[i])
    eval(parse(text=
                 paste0("attributes(", iPaste, ") <- att")
    ))

    # Add leaves
    leaves <- matCollapse[matCollapse == i]
    if(length(leaves) > 0){
      for(l in 1:length(leaves)){

        # Add element
        lPaste <- paste0(iPaste, "[[", l, "]]")
        eval(parse(text=
                     paste0(lPaste, "<-list()")
        ))

        # Add attributes
        members <- 1
        leaf <- names(leaves)[l]
        height <- 0
        att <- list(members = members, height = height, label = leaf, leaf=TRUE)
        eval(parse(text=
                     paste0("attributes(", lPaste, ") <- att")
        ))
      }
    }
  }

  class(aList) <- "dendrogram"
  return(aList)
}

doSplit <- function(Dfram, nBoots = 200, nGenes = nrow(Dfram)/10, seed = 123, agg = "ward.D2"){

  # Order the chemicals in alphbetical order to generate predictable splits
  Dfram <- Dfram[,order(colnames(Dfram))]

  # Get consensus module
  set.seed(123)
  modVec <- do.call(c, lapply(1:nBoots, function(x){
    Dboot <- Dfram[sample(nrow(Dfram), nrow(Dfram), replace = TRUE),]
    MEANs <- sort(apply(Dboot, 1, function(x) sum(x^2)), decreasing = T)
    MEANs <- MEANs[1:(nGenes)]
    Dboot <- Dboot[names(MEANs),]
    dDist <- dist(t(Dboot), method = "euclidean")
    dClust <- hcopt(dDist, method = agg)
    mods <- paste(cutree(dClust, k = 2), collapse = "")
    return(mods)
  }))
  modTab <- sort(table(modVec), decreasing = T)
  modList <- strsplit(names(modTab), "")

  mods <- as.factor(modList[[1]]); names(mods) <- colnames(Dfram)
  modProp <- modTab[1]/nBoots; names(modProp) <- NULL

  # Get stability statistics for each sample
  modFram <- do.call(rbind, lapply(modList, function(x) x)); colnames(modFram) <- colnames(Dfram)

  modStab <- do.call(c, lapply(1:2, function(x){
    mSub <- modFram[,mods == x, drop = FALSE]
    mat <- matrix(NA, ncol(mSub), ncol(mSub))
    for(i in 1:ncol(mSub)) for(j in 1:ncol(mSub)) mat[i,j] <- sum(as.numeric(mSub[,i] == mSub[,j])*modTab/nBoots)
    diag(mat) <- NA
    out <- colMeans(mat, na.rm = T)
    names(out) <- colnames(mSub)
    out
  }))[colnames(Dfram)]

  if(min(table(mods)) == 1) mods <- NULL

  return(list(mods = mods, propBoots = modProp, clustStab = modStab))
}


runHClust <- function(Dfram, nGenes = nrow(Dfram)/10){
  # Filter by signficance
  MEANs <- sort(apply(Dfram, 1, function(x) sum(x^2)), decreasing = T)
  MEANs <- MEANs[1:(nGenes)]
  Dfram <- Dfram[names(MEANs),order(colnames(Dfram))]
  dDist <- dist(t(Dfram), method = "euclidean")
  dClust <- hcopt(dDist, method = "ward.D2")
  return(dClust)
}

runK2Clust <- function(Dfram, nGenes = nrow(Dfram)/10, agg = "ward.D2"){

  # Create Splits of the data
  taxList <- list(list(colnames(Dfram)))
  resList <- stabList <- list(list(NULL))
  iter <- 1
  while(max(unlist(lapply(taxList[[iter]], function(x) length(x[x!="VEHICLE"])))) > 2){
    outList <- lapply(taxList[[iter]], function(samps){
      if(length(samps)>2){
        Dsub <- Dfram[,samps]
        outList <- doSplit(Dsub, nGenes = nGenes, agg = agg)
      } else {
        outList <- list(mods = samps, propBoots = NULL, clustStab = NULL)
      }
      return(outList)
    })
    iter <- iter+1
    taxList[[iter]] <- list()
    outMods <- lapply(outList, function(x) x[[1]])
    slot <- 1
    for(i in 1:length(outMods)){
      taxList[[iter]][[slot]] <- names(outMods[[i]])[outMods[[i]] == 1]; slot <- slot + 1
      taxList[[iter]][[slot]] <- names(outMods[[i]])[outMods[[i]] == 2]; slot <- slot + 1
    }

    if(length(taxList[[iter]])==0){
      taxList[[iter]] <- list(NULL)
    }

    resList[[iter]] <- list()
    stabList[[iter]] <- list()
    outRes <- lapply(outList, function(x) x[[2]])
    outStab <- lapply(outList, function(x) x[[3]])
    slot <- 1
    for(i in 1:length(outRes)){
      resList[[iter]][[slot]] <- outRes[[i]]
      stabList[[iter]][[slot]] <- outStab[[i]]
      slot <- slot + 1
    }
  }
  resList <- resList[-1]
  stabList <- stabList[-1]

  if(is.null(unlist(taxList[[iter]]))) taxList <- taxList[-iter]

  # Get instances where the split had > 2 samples in cluster
  modList <- lapply(taxList[-1], function(x){
    combs <- 1:(length(x)/2)
    combList <- list()
    for(i in combs){
      combList[[i]] <- list()
      combList[[i]][[1]] <- x[[i*2-1]]
      combList[[i]][[2]] <- x[[i*2]]
      if(length(combList[[i]]) > 0){
        if(length(combList[[i]][[1]]) < 2 | length(combList[[i]][[2]]) < 2) combList[[i]] <- list()
      }
    }
    combList
  })

  # Create list of modules and bootstrap stats
  finList <- lapply(1:length(modList), function(x){
    lapply(1:length(modList[[x]]), function(y){
      modSub <- modList[[x]][[y]]
      if(length(modSub) > 0){
        resSub <- resList[[x]][[y]]
        stabSub <- stabList[[x]][[y]]
        return(list(chems = modSub, bootP = resSub, sampStab = stabSub))
      } else {
        return(NULL)
      }
    })
  })
  finUnList <- unlist(finList, recursive = F)
  finUnList <- finUnList[!unlist(lapply(finUnList, is.null))]

  return(finUnList)

}

# Run DGE on clusters
runDGE_mods <- function(eSet, finUnList, Dfram, nGenes = nrow(Dfram)/10){
  finUnList <- lapply(finUnList, function(x){

    # Filter by signficance
    Dsub <- Dfram[,c(x[[1]][[1]], x[[1]][[2]])]
    MEANs <- sort(apply(Dsub, 1, function(x) sum(x^2)), decreasing = T)
    MEANs <- MEANs[1:nGenes]

    # Run DGE

    # Create module variable
    mods <- as.factor(c(rep(1, length(x[[1]][[1]])), rep(2, length(x[[1]][[2]]))))
    names(mods) <- paste0("X", c(x[[1]][[1]], x[[1]][[2]]))
    mods <- mods[names(mods) != "XVEHICLE"]
    eSet$Chemical <- paste0("X", eSet$Chemical)
    eSetSub <- eSet[names(MEANs), eSet$Chemical %in% names(mods)]
    eSetSub$mod <- mods[eSetSub$Chemical]

    # Create design matrix
    design <- model.matrix(~ 0 + Chemical, data = pData(eSetSub))
    colnames(design) <- sub("Chemical", "", colnames(design))

    # Create contrast matrix
    dLeft <- paste0("(",paste(names(mods)[mods == 2], collapse = "+"), ")/", sum(mods == 2))
    dRight <- paste0("(",paste(names(mods)[mods == 1], collapse = "+"), ")/", sum(mods == 1))
    dFull <- paste(dLeft, dRight, sep = "-")
    Contrast <- makeContrasts(contrasts=dFull, levels = colnames(design))

    print(mods)

    # Create contrasts between chemicals
    corfit <- duplicateCorrelation(eSetSub, design, block = eSetSub$Plate)

    if(corfit$consensus.correlation %in% c(-1, 1) | is.na(corfit$consensus.correlation)){
      suppressWarnings(fit <- lmFit(eSetSub, design))
    } else {
      suppressWarnings(fit <- lmFit(eSetSub, design, correlation = corfit$consensus.correlation, block = eSetSub$Plate))
    }
    suppressWarnings(fit <- contrasts.fit(fit, contrasts = Contrast))
    suppressWarnings(fit <- eBayes(fit, trend=TRUE))

    # Add DGE Results
    x$res <- topTable(fit, number = Inf, sort.by = "P")

    # Add mean z-scores to results
    x$res$mean_z1 <- rowMeans(Dfram[rownames(x$res),x$chems[[1]][x$chems[[1]]!="VEHICLE"]])
    x$res$mean_z2 <- rowMeans(Dfram[rownames(x$res),x$chems[[2]][x$chems[[2]]!="VEHICLE"]])
    return(x)
  })
}

# Plot Dendrogram

plotDendro <- function(dClust, splitNames = LETTERS[1:length(finUnList)]){
  dend_segs <- dendro_data(dClust, type = "rectangle")$segments
  dend_horiz <- unique(dend_segs[dend_segs$y == dend_segs$yend & dend_segs$y > 1, c("x", "y")])
  dend_horiz <- dend_horiz[order(dend_horiz$x),]
  dend_horiz <- dend_horiz[order(dend_horiz$y, decreasing = T),]
  dend_horiz$label <- splitNames

  # Add Labels of leaves
  leaf_labels = get_leaves_attr(dClust, "label")

  # Create plot of dendrogram
  dendroP <- ggdendrogram(dClust) +
    geom_label(data = dend_horiz, aes(x = x, y = y, label = label), size = 6, color = "red") +
    theme_bw() +
    scale_x_reverse(breaks = 1:length(leaf_labels), label = leaf_labels, position = "top") +
    scale_y_reverse() +
    coord_flip() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(hjust = 0, size = 18),
      axis.text.x = element_blank()
    )
  return(dendroP)
}

# Fix FDRs
fixFDRs <- function(finUnList){
    dgeFram <- do.call(rbind, lapply(1:length(finUnList), function(x){
        res <- finUnList[[x]]$res
        res$slot <- x
        return(res)
      }))
      dgeFram$adj.P.Val <- p.adjust(dgeFram$P.Value, method = "BH")
      outList <- lapply(1:length(finUnList), function(x){
        xObj <- finUnList[[x]]
        RN <- rownames(xObj$res)
        xObj$res <- dgeFram[dgeFram$slot == x,colnames(xObj$res)]

        rownames(xObj$res) <- RN
        xObj
      })
  return(outList)
}

# Generate network diagram
generateVisNetwork <- function(finUnList){

  require(visNetwork)

  # Generate Matrix from visNetwork
  mat <- matrix(0, nrow = length(finUnList), ncol = ncol(Dfram))
  colnames(mat) <- colnames(Dfram)
  for(i in 1:length(finUnList)){
    mat[i,finUnList[[i]]$chems[[1]]] <- 1
    mat[i,finUnList[[i]]$chems[[2]]] <- 2
  }

  rownames(mat) <- LETTERS[1:nrow(mat)]

  # Calculate sizes
  sizes <- apply(mat, 1, function(x) sum(x!=0))

  # Add Labels
  titles <- apply(mat, 1, function(x) paste(colnames(mat)[x!=0], collapse = "<br>")); names(titles) <- names(sizes)

  source <- c()
  target <- c()
  k <- 1

  # Get edges
  for(i in 1:nrow(mat)){

    source <- c(source, rep(rownames(mat)[i], 2))
    matRow <- mat[i,]
    sub1 <- which(matRow == 1)[1]
    sub2 <- which(matRow == 2)[1]
    matSub1 <- mat[-(1:i), sub1];names(matSub1) <- rownames(mat)[-(1:i)]
    matSub2 <- mat[-(1:i), sub2];names(matSub2) <- rownames(mat)[-(1:i)]

    target1 <- names(matSub1)[which(matSub1!=0)[1]]
    if(is.na(target1)){
      target1 <- letters[k]
      sizes[target1] <- sum(matRow == 1)
      titles[target1] <- paste(colnames(mat)[matRow == 1], collapse = "<br>")
      k <- k+1}
    target2 <- names(matSub1)[which(matSub2!=0)[1]]
    if(is.na(target2)){
      target2 <- letters[k]
      sizes[target2] <- sum(matRow == 2)
      titles[target2] <- paste(colnames(mat)[matRow == 2], collapse = "<br>")
      k <- k+1}
    target <- c(target, target1, target2)
  }

  # Set terminal nodes to 0 sizes
  sizes[names(sizes) %in% letters] <- 0

  # Add Labels
  labs <- titles; labs <- gsub("<br>", "\n", labs); labs[names(labs) %in% LETTERS] <- names(labs)[names(labs) %in% LETTERS]

  # Add shapes
  shapes <- rep("diamond", length(sizes)); names(shapes) <- names(sizes); shapes[names(shapes) %in% letters] <- "box"

  # Get terminal node level
  levs <- sapply(2:nrow(mat), function(x){
    matSub <- mat[1:x,]
    matSub <- matSub[,matSub[x,]!=0]
    matSub <- matSub[,1]
    sum(matSub!=0)
  })
  levs <- c(1, levs)
  levs <- c(levs, rep(max(levs) + 1, length(sizes) - length(levs)))
  names(levs) <- names(sizes)

  # Genreate network
  nodes <- data.frame(id = names(levs), title = titles, level = levs, label = labs, shape = shapes)
  edges <- data.frame(from = source, to = target)

  # Generate plot
  p <- visNetwork(nodes, edges) %>%
    visEdges(arrows = "to") %>%
    visHierarchicalLayout(direction = "LR")

  return(p)
}
