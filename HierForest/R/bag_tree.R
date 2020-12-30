bag_tree <- function(features, group, class, featPtree){
  
  # Bag samples and randomly choose features
  bootSamps <- sample(nrow(features), replace = T)
  randFeats <- sample(ncol(features), featPtree, replace = F)

  # Subset features groups and classes
  fSamp <- features[bootSamps, randFeats]
  gSamp <- group[bootSamps]
  cSamp <- class[bootSamps]
  
  # Report samples used in bootstrap
  samp_included = rownames(features)[unique(bootSamps)]

  return(list(
    features = fSamp,
    group = gSamp,
    class = cSamp,
    samp_included = samp_included
  ))
}
