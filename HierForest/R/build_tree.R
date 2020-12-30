build_tree <- function(features, class, max.depth = Inf, min.nodesize = 1, samp_included = NULL){
  
  
  # Initialize
  tree <- list()
  tree[[1]] <- wrap_split(features, class)
  tree[[1]]$dirs <- 1
  tree[[1]]$fUsed <- names(tree[[1]]$gini)
  tree[[1]]$samp_included <- samp_included
  treeOut <- tree
  mDepth <- 1
  
  # Make sure tree isn't built beyond number of features
  if(max.depth == Inf) max.depth <- ncol(features) - 1
  
  # Build tree to max depth
  while(mDepth <= max.depth & length(unlist(treeOut)) > 0){
    outList <- lapply(treeOut, function(split){
      if(!is.null(split) & length(split$group1Class) > min.nodesize & sum(split$group1Class != split$group1Consensus) > 0){
        branch1 <- wrap_branch(split$group1, features, class, split$fUsed)
        branch1$dirs <- c(split$dirs, 1)
        branch1$fUsed <- c(split$fUsed, names(branch1$gini))
      } else {
        branch1 <- NULL
      }
      if(!is.null(split) & length(split$group2Class) > min.nodesize & sum(split$group2Class != split$group2Consensus) > 0){
        branch2 <- wrap_branch(split$group2, features, class, split$fUsed)
        branch2$dirs <- c(split$dirs, 2)
        branch2$fUsed <- c(split$fUsed, names(branch2$gini))
      } else {
        branch2 <- NULL
      }
      outBranch <- list(branch1, branch2)
      outBranch <- outBranch[!unlist(lapply(outBranch, is.null))]
      if(length(outBranch) > 0) return(outBranch) else return(NULL)
    })
    treeOut <- unlist(outList, recursive = F)
    mDepth <- mDepth + 1
    tree[[mDepth]] <- treeOut
  }
  return(tree)
}
