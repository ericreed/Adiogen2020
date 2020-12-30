wrap_split <- function(features, class){

  # Calculate best split
  gSplit <- split_groups(features, class)

  # Return best split
  return(gSplit)

}
