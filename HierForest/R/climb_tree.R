climb_tree <- function(samp, tree){
  branch <- 1
  direct <- 1
  depth <- length(tree)
  whichNode <- 999
  while(branch <= depth & length(whichNode) != 0){
    split <- tree[[branch]]
    if(branch == 1){
      node <- split
      dirI <- as.numeric(samp[names(node$gini)] > node$value) + 1
      direct <- c(direct, dirI)
      branch <- branch + 1
    } else {
      whichNode <- which(
        unlist(
          lapply(
            lapply(split, function(x) x$dirs), function(x) identical(x, direct)
          )
        )
      )
      if(length(whichNode) != 0){
        node <- split[[whichNode]]
        dirI <- as.numeric(samp[names(node$gini)] > node$value) + 1
        direct <- c(direct, dirI)
        branch <- branch + 1
      }
    }
  }
  termDir <- paste0("group", direct[length(direct)], "Consensus")
  class <- node[[termDir]]
  return(class)
}
