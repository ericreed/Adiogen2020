require(Biobase)

# Directories
datDir <- "../data"
resDir <- "../results"
supDir <- "support"

## Subdirectories
esetDir <- file.path(datDir, "eSets")
cvDataDir <- file.path(datDir, "cvDat")

# Files to read in
eSetFile <- file.path(esetDir, paste0("Obesogen_3DGE_ExpA_AllSamples_Combat_LERem_", "011a", ".rds"))

# Files to write out
outDat <- file.path(cvDataDir, "repCVdata.rds")

# Read in helper functions
source(file.path(supDir, "preproc_functions.R"))

# Read in eSet
eSet <- readRDS(eSetFile)

# Format for random forest
labList <- format4rn(eSet)

# Create 10 fold CV group list (20 iterations)
kCV <- 10
nReps <- 10

groups <- unique(labList$group)
index <- rep(1:kCV, length.out = length(groups))
set.seed(123)
iterList <- lapply(1:nReps, function(i){
    index <- sample(index)
    lapply(1:kCV, function(x) groups[index == x])
})

# Combine into list
outList <- list(labList = labList, iterList = iterList)

# Save
saveRDS(outList, outDat)
