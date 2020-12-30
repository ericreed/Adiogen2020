#####MOUSE-MEF Phosphorylation Data Set#####
##Affymetrix Mouse Genome 430A 2.0 Array
source("parameters.R")

library(GEOquery)
library(ArrayExpress)
library(frma)
library(biomaRt)
library(reshape2)
library(limma)
library(simpleaffy)
library(RColorBrewer)
library(affyPLM)
library(rgl)
library(annotate)
library(CBMRtools)
library(gplots)
library(oligo)
library(pdInfoBuilder)
library(ff)
library(doMC)
registerDoMC(1)
require(XML)
require(hgu219.db)

# Create annotation
download.file("http://mbni.org/customcdf/22.0.0/ensg.download/HGU219_Hs_ENSG_22.0.0.zip",file.path(HMAdir, "tmp.zip"))
unzip(file.path(HMAdir, "tmp.zip"), exdir = HMAdir)
file.remove(file.path(HMAdir, "tmp.zip"))

z <- cdf2table(file.path(HMAdir, "HGU219_Hs_ENSG.cdf"))

seed <- new("GenericPDInfoPkgSeed",
             table=z,
             author = "Eric Reed",
             email = "reeder@bu.edu",
             species = "homo sapients",
             pkgName = "pd.HGU219.Hs.ENSG")

makePdInfoPackage(seed, HMAdir)
install.packages(file.path(HMAdir, "pd.HGU219.Hs.ENSG"), repos = NULL, type = "source")

# Get list of cel files
cels <- list.files(file.path(HMAdir, "GSE70353_RAW"), pattern = "cel")

#loading and normalizing data
RawAffy <- read.celfiles(file.path(HMAdir, "GSE70353_RAW", cels), pkgname = "pd.HGU219.Hs.ENSG")

#perform RMA
eSet = rma(RawAffy)
colnames(eSet) <- substr(colnames(eSet), 1, 10)

# Create phenotype data

### Get sample information from cel file names
iLines <- readLines(file.path(HMAdir, "GSE70353_series_matrix.txt"))
iLines <- iLines[grepl("Sample_geo_accession|Sample_characteristics", iLines)]
iLines <- gsub("!Sample_geo_accession|!Sample_characteristics_ch1|\t", "", iLines)
iLines <- gsub("\"\"", "\"", iLines)
iLines <- gsub("\"", ":", iLines)
iLines <- strsplit(iLines, ":")
iLines <- lapply(iLines, function(x) x[-1])

### Format and get column name
IDs <- iLines[[1]]; iLines <- iLines[-1]
varNames <- unlist(lapply(iLines, function(x) x[1]))
valList <- lapply(iLines, function(x) trimws(x[!x %in% varNames]))
varNames <- trimws(varNames)

### Create data frame of pData
pdata <- do.call(cbind, valList); colnames(pdata) <- varNames
pdata <- data.frame(ID = IDs, pdata, row.names = IDs, stringsAsFactors = F)
pdata <- pdata[colnames(eSet),]

# Format columns to numeric and factors
numCols <- c(4:27)
for(i in numCols) pdata[,i] <- as.numeric(pdata[,i])

# Create feature data set
fDat <- fData(eSet)
fDat$ENSEMBL <- sub("_at", "", rownames(fDat))

# Add symbols
toHGNC <- select(hgu219.db,
                 key=fDat$ENSEMBL,
                 columns=c("ENSEMBL", "SYMBOL"),
                 keytype="ENSEMBL")
toHGNC <- toHGNC[!duplicated(toHGNC$ENSEMBL),]
rownames(toHGNC) <- toHGNC$ENSEMBL
toHGNC <- toHGNC[fDat$ENSEMBL,]
if(identical(rownames(toHGNC), fDat$ENSEMBL)) fDat$SYMBOL <- toHGNC$SYMBOL

# Make ESET
eSet = to.eSet(exprs(eSet), pdata, fDat)

# Make symbols the rownames
eSet <- eSet[!duplicated(fData(eSet)$SYMBOL),]
eSet <- eSet[!is.na(fData(eSet)$SYMBOL),]
rownames(eSet) <- fData(eSet)$SYMBOL

# Save
saveRDS(eSet, eSet_Human_124)
