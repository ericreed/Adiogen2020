
GiniImportance <- function(forest){

  giniVec <- unlist( lapply(forest, function(tree){

    # Make the structure consistent
    tree[[1]][[1]] <- tree[[1]]
    tree[[1]][2:length(tree[[1]])] <- NULL; names(tree[[1]]) <- NULL

    # Get weighted gini importance
    giniWeighted <- unlist(lapply(tree, function(lev){
      giniSplit <- lapply(lev, function(s){
        sGini <- (1-s$gini) * sum(length(s$group1), length(s$group2))
      })
    }))

    # Sum weights of features used twice (possible for divergent samples)
    namDup <- unique(names(giniWeighted)[duplicated(names(giniWeighted))])
    if(length(namDup)>0){
      gRep <- sapply(namDup, function(x) sum(giniWeighted[names(giniWeighted)==x]))
      giniWeighted <- giniWeighted[!names(giniWeighted) %in% namDup]
      giniWeighted <- c(giniWeighted, gRep)
    }

    return(giniWeighted)
  }))

  # Get mean of gini importance
  geneImp <- sort(sapply(unique(names(giniVec)), function(x) sum(giniVec[names(giniVec) == x])), decreasing = T)/length(forest)

  return(geneImp)

}
