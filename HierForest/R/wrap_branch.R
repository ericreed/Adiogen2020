wrap_branch <- function(groups, features, class, fUsed){

  # Subset for branch to run split
  gKeep <- rownames(features) %in% groups
  fSub <- features[gKeep,!colnames(features) %in% fUsed, drop = FALSE]
  cSub <- class[gKeep]

  return(wrap_split(fSub, cSub))

}
