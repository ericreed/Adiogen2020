get_votes <- function(tree, gMeans){
  votes <- apply(gMeans, 1, climb_tree, tree)
  return(votes)
}
