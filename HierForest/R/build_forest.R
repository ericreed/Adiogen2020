build_forest <- function(features, group, class, max.depth = Inf, min.nodesize = 1, featPtree = 100, nTrees = 2000, cores = 1, mergeFunction = groupMeans, bagFunction = bag_tree){

  #max.depth = Inf; featPtree = sqrt(ncol(features)); nTrees = 10; cores = 1; mergeFunction = NoTransform; bagFunction = bag_tree;min.nodesize = 1
  
  # Repeat tree bagging
  treeList <- mclapply(1:nTrees, function(i){
    bagList <- bagFunction(features, group, class, featPtree)
    datMerge <- mergeFunction(bagList$features, bagList$group, bagList$class)
    build_tree(datMerge$features, datMerge$class, max.depth, min.nodesize, bagList$samp_included)
  }, mc.cores = cores, mc.preschedule=F)

  return(treeList)

}
