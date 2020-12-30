oob_votes <- function(tree, features, group, mergeFunction){
  
  # Get samples not in tree
  samps_use <- rownames(features)[!rownames(features) %in% tree[[1]]$samp_included]
  
  # Subset features, class, and group
  keep <- rownames(features) %in% samps_use
  fSub <- as.matrix(features[keep,])
  gSub <- group[keep]
  
  # Transform data
  gMeans <- mergeFunction(fSub, gSub, class = NA)
  
  # Get votes
  votes <- apply(gMeans$features, 1, climb_tree, tree)
 
  # Return votes
  return(votes)
}
