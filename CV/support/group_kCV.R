# Run Cross Validation
group_kCV <- function(Flist,
                 iterList,
                 preproc_train_function = NULL,
                 preproc_test_function = NULL,
                 bag_function,
                 transform_function = NULL,
                 split_metric,
                 nTrees = 2000,
                 nCores = 10,
                 max.depth = Inf,
                 file_prefix
){

  reps <- length(iterList)
  outRepeated <- lapply(1:reps, function(i, file_prefix){
    diditwork <- try({
      print(paste0("Iteration = ", i))
      groupList <- iterList[[i]]
      kCVout <- lapply(groupList, function(groupVec){

        # Remove group k
        keep <- !Flist$group %in% groupVec
        train <- lapply(Flist, function(x){
          if(!is.null(dim(x))) x <- x[keep,] else x <- x[keep]
          return(x)
        })

        # Preproccess training data
        if(!is.null(preproc_train_function)) train <- preproc_train_function(train)

        # Train Forest
        forest_train <- build_forest(train$features,
                                     train$groups,
                                     train$class,
                                     max.depth = max.depth,
                                     featPtree = sqrt(ncol(train$features)),
                                     nTrees = nTrees,
                                     cores = nCores,
                                     bagFunction = bag_function,
                                     mergeFunction = transform_function)

        # Get OOB Estimates
        oobEstimates <- oob_estimate(forest_train, train$features, train$class, train$groups, mergeFunction = transform_function)

        if(is.na(oobEstimates$class[1])) {
          oobEstimates$class <- train$class[rownames(oobEstimates)]
          oobEstimates$group <- train$groups[rownames(oobEstimates)]
        }

        # Test Forest
        test <- lapply(Flist, function(x){
          if(!is.null(dim(x))) x <- x[!keep,] else x <- x[!keep]
          return(x)
        })

        # Preproccess test data
        if(!is.null(preproc_test_function)) test <- preproc_test_function(test)

        # Predict test data
        testEstimates <- predict_forest(test$features, test$groups, forest_train, mergeFunction = transform_function)

        return(list(
          testEst = testEstimates,
          oobVoting = oobEstimates)
        )
      })

      saveRDS(kCVout, paste0(file_prefix, gsub("-| |:", "_", Sys.time()), "_", i, ".rds"))
    }, silent = TRUE)
    if(class(diditwork) != "try-error") return("Success") else return("Failed")
  }, file_prefix = file_prefix)
  return(outRepeated)
}
