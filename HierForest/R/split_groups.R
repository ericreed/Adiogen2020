# Define split for each row
split_groups <- function(features, class){

  # Get best split for every gene
  splitMat <- matrix(1, length(class), length(class)); splitMat[lower.tri(splitMat)] <- 0
  splitStats <- apply(features, 2, function(feat){
    ord <- order(feat)
    feat <- feat[ord]; class <- class[ord]
    gini_index(class, splitMat = splitMat, rMin = T)
  }); names(splitStats) <- colnames(features); splitStats <- sort(splitStats)
  
  # Get split for best gene
  gBestSplit <- splitStats[1]
  feat <- features[, names(gBestSplit)]
  ord <- order(feat)
  feat <- feat[ord]; class <- class[ord]
  gini_vec <- gini_index(class, splitMat = splitMat)
  splitVal <- feat[which.min(gini_vec)]

  return(list(value = splitVal,
              gini = gBestSplit,
              group1 = names(feat)[feat<=splitVal],
              group1Class = class[feat<=splitVal],
              group1Consensus = names(sort(table(class[feat<=splitVal])))[2],
              group2 = names(feat)[feat>splitVal],
              group2Class = class[feat>splitVal],
              group2Consensus = names(sort(table(class[feat>splitVal])))[2]
  ))
}
