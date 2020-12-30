oob_estimate <- function(forest, features, class, group, mergeFunction = groupMeans){
  
  # Get predictions
  predList <- lapply(forest, oob_votes, features, group, mergeFunction = mergeFunction)
  
  # Now build matrix of results
  gNames <- sort(unique(unlist(lapply(predList, names))))
  predMat <- matrix(NA, nrow = length(predList), ncol = length(gNames)); colnames(predMat) <- gNames
  for(i in 1:length(predList)) predMat[i,names(predList[[i]])] <- predList[[i]]
  
  # Covert to data.frame and get colSum
  predFram <- as.data.frame(predMat)
  results <- colMeans(predFram == "Yes", na.rm = T)
  
  # Get classes
  gUnique <- unique(data.frame(class = class, group = group)); rownames(gUnique) <- gUnique$group
  
  # Return data fram of OOB results and estimates
  resFram <- data.frame(votes = results, class = gUnique[names(results), "class"], row.names = names(results))
  
  return(resFram)

}
