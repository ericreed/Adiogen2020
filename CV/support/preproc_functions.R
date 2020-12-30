# Format expression set for random forest CV
format4rn <- function(eSet) {

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
  
  return(labList)
}


# Adip data preproc function (Bag Merge)
preproc_train_adip <- function(Flist){
  
  #Filter Genes
  runANOVA <- function(gene){
    df <- data.frame(e = Flist$features[,gene], g = Flist$group)
    m <- anova(lm(e~g, data = df))$`Pr(>F)`[1]
    return(m)
  }
  
  Aout <- do.call(c, lapply(colnames(Flist$features), runANOVA))
  names(Aout) <- colnames(Flist$features)
  Aout <- sort(Aout); Aout <- p.adjust(Aout, method = "BH")
  Flist$features <- Flist$features[,names(Aout)[Aout < 0.2]]
  
  return(Flist)
}

# Adip data preproc function (Pre-merge)
preproc_train_adip_premerge <- function(Flist){
  
  # Filter Genes
  runANOVA <- function(gene){
    df <- data.frame(e = Flist$features[,gene], g = Flist$group)
    m <- anova(lm(e~g, data = df))$`Pr(>F)`[1]
    return(m)
  }
  
  Aout <- do.call(c, lapply(colnames(Flist$features), runANOVA))
  names(Aout) <- colnames(Flist$features)
  Aout <- sort(Aout); Aout <- p.adjust(Aout, method = "BH")
  Flist$features <- Flist$features[,names(Aout)[Aout < 0.2]]
  
  # Tranform genes
  Flist <- groupMeans(Flist$features, Flist$groups, Flist$class)
  Flist$groups <- rownames(Flist$features)
  
  return(Flist)
}

# Adip data preproc function (Pre-merge)
preproc_test_adip_premerge <- function(Flist){
  
  # Tranform genes
  Flist <- groupMeans(Flist$features, Flist$groups, Flist$class)
  Flist$groups <- rownames(Flist$features)
  
  return(Flist)
}