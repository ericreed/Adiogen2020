require(parallel)

# Directories
datDir <- "../data"
resDir <- "../results"
supDir <- "support"

## Subdirectories
esetDir <- file.path(datDir, "eSets")
outputDir <- file.path(resDir, "CVoutput")
cvDataDir <- file.path(datDir, "cvDat")

packDir <- "../HierForest/"

# Files to read in
cvDatFile <- file.path(cvDataDir, "repCVdata.rds")

# Files to write out
outFile <- file.path(outputDir, "AdipMean_repCV_001a_BagMerge_")

# Read in package functions
thisDir <- getwd()
setwd(packDir)
source("R/loadFunctions.R")
setwd(thisDir)

# Read in helper functions
source(file.path(supDir, "preproc_functions.R"))
source(file.path(supDir, "group_kCV.R"))

# Format for random forest
cvDat <- readRDS(cvDatFile)
labList <- cvDat$labList
iterList <- cvDat$iterList

preproc_train_function = preproc_train_adip
preproc_test_function = NULL
bag_function = bag_tree
transform_function = groupMeans
split_metric  = gini_index
nTrees = 2000
nCores = 15
max.depth = Inf

# Run cv
group_kCV(labList,
          iterList,
          preproc_train_function = preproc_train_function,
          preproc_test_function = preproc_test_function,
          bag_function = bag_function,
          transform_function = transform_function,
          split_metric = split_metric,
          nTrees = nTrees,
          nCores = nCores,
          max.depth = max.depth,
          file_prefix = outFile
)
