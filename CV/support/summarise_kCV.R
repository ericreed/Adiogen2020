# Summarise repeated k-fold-cross validation
summarise_kCV <- function(kCVout,
                          threshFunction = oob_thresh_f1,
                          fTrees = seq(0.05, 1, by = 0.05),
                          collapse = F){

  kCV_list <- lapply(kCVout,
                     summarise_kCVi,
                     fTrees = fTrees,
                     threshFunction = threshFunction,
                     collapse = collapse)

  return(kCV_list)

}


# Summarise a single k-fold cross-valudation
summarise_kCVi <- function(kCVi,
                           threshFunction = oob_thresh_f1,
                           collapse = F
                           ){

  # Create list object to output
  CVsum <- list()

  # Add get true class and groups
  classes <- unique(do.call(rbind, lapply(kCVi, function(x) {
    xDF <- x$oobVoting[, colnames(x$oobVoting) != "votes", drop = FALSE]
    xDF$Obs <- rownames(xDF)
    rownames(xDF) <- NULL
    return(xDF)
  }))); rownames(classes) <- classes$Obs

  # Collapse
  if(collapse) {
    kCVi <- lapply(kCVi, Collapse, classes)
    classes <- unique(classes[, c("class", "group")])
    rownames(classes) <- classes$group
  }

  # Get oob thresholds to test estimates
  kCVi <- lapply(kCVi, threshFunction)

  # Get prediction frame
  predFram <- do.call(rbind, lapply(kCVi, function(x) {
    data.frame(votes = x$testEst, oob_thresh = x$oob_thresh, row.names = names(x$testEst))
  }))
  predFram <- predFram[order(predFram$votes),]

  # Add prediction
  predFram$pred <- c("No", "Yes")[as.numeric(predFram$votes > predFram$oob_thresh) + 1]

  predFram$class <- classes[rownames(predFram), "class"]

  CVsum$predFram <- predFram

  # Add performance
  CVsum$perfFram <- getPerfStats(predFram)

  return(CVsum)
}

oob_thresh_f1 <- function(kCVii, collapse = F){

  oob_res <- kCVii$oobVoting[order(kCVii$oobVoting$votes),]

  oob_res$F1 <- sapply(1:nrow(oob_res), function(x){

    # Caulcate f1
    sens <- sum(oob_res$class[x:nrow(oob_res)] == "Yes") / sum(oob_res$class == "Yes")
    prec <- sum(oob_res$class[x:nrow(oob_res)] == "Yes") / length(x:nrow(oob_res))
    f1 <- sum(na.omit(c(2 * ((sens*prec)/(sens + prec)), 0)), na.rm = F)

    return(f1)
  })

  #kCVii$oob_thresh <- oob_res$votes[which.max(oob_res$F1) - 1]
  kCVii$oob_thresh <- mean(oob_res$votes[c(which.max(oob_res$F1), which.max(oob_res$F1) - 1)])

  return(kCVii)
}

getPerfStats <- function(predFram){

  # Get AUC
  auc <- roc(predFram$class, predFram$votes)$auc

  # Get Sens, Spec, Prec, Acc, F1, BalAcc
  sens <- sum(predFram$class == "Yes" & predFram$pred == "Yes") / sum(predFram$class == "Yes")
  spec <- sum(predFram$class == "No" & predFram$pred == "No") / sum(predFram$class == "No")
  prec <- sum(predFram$class == "Yes" & predFram$pred == "Yes") / sum(predFram$pred == "Yes")
  acc <- mean(predFram$class == predFram$pred)
  f1 <- 2 * (sens*prec)/(sens + prec)
  balacc <- mean(c(sens, spec))

  # Format as data.frame (Easy to combine to other results)
  perfFram <- data.frame(auc, balacc, f1, prec, sens, spec, acc)

  return(perfFram)
}

# Function to collapse group voting
Collapse <- function(kCVii, classes){

  require(dplyr)

  # Get test estimates
  testEst <- kCVii$testEst

  # Create data.frame of test estimates
  testEst <- data.frame(votes = kCVii$testEst,
                        group = classes[names(kCVii$testEst), "group"])

  # Collapse test estimates
  testEst <- as.data.frame(
    testEst %>%
      group_by(group) %>%
      summarise(votes = mean(votes))
  )
  testNames <- as.character(testEst$group)
  testEst <- testEst$votes
  names(testEst) <- testNames
  kCVii$testEst <- testEst


  # Collapse oobVoting
  oobVoting <- kCVii$oobVoting
  oobVoting <- as.data.frame(
    oobVoting %>%
      group_by(group, class) %>%
      summarise(votes = mean(votes))
  )
  rownames(oobVoting) <- as.character(oobVoting$group)
  kCVii$oobVoting <- oobVoting[, c("votes", "class")]

  return(kCVii)
}
