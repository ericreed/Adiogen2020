predict_forest <- function(feature_test, group_test, forest, mergeFunction){

  gMeans <- mergeFunction(feature_test, group_test, class = NA)
  voteFrame <- as.data.frame(do.call(rbind, lapply(forest, get_votes, gMeans$features)))
  results <- colMeans(voteFrame == "Yes")
  return(results)

}
