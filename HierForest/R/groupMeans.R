groupMeans <- function(features, group, class){

  # Create data matrix of just group info
  gInfo <- data.frame(group, class)
  gUnique <- unique(gInfo); rownames(gUnique) <- gUnique$group

  group <- as.factor(as.character(group))

  # Create X matrix
  X <- model.matrix(~ 0 + group)
  colnames(X) <- sub("group", "", colnames(X))

  # Calculate means for each chemical
  Gmeans <- solve( t(X) %*% X ) %*% t(X) %*% features

  # Sort gUnique
  gUnique <- gUnique[rownames(Gmeans),]

  # Return List of Features and Classes
  return(list(features = Gmeans, class = gUnique$class))

}
