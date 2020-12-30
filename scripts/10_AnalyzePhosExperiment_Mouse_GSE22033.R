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
library(gplots)
library(oligo)
library(pdInfoBuilder)
library(ff)
library(doMC)
registerDoMC(1)
require(XML)
require(mouse430a2.db)

# Create annotation
# download.file("http://mbni.org/customcdf/22.0.0/ensg.download/Mouse430A2_Mm_ENSG_22.0.0.zip",file.path(MPdir, "tmp.zip"))
# unzip(file.path(MPdir, "tmp.zip"), exdir = MPdir)
# file.remove(file.path(MPdir, "tmp.zip"))
# 
# z <- cdf2table(file.path(MPdir, "Mouse430A2_Mm_ENSG.cdf"))
# 
# seed <- new("GenericPDInfoPkgSeed",
#             table=z,
#             author = "Eric Reed",
#             email = "reeder@bu.edu",
#             species = "mus musculus",
#             pkgName = "pd.Mouse430A2.Mm.ENSG")
# 
# makePdInfoPackage(seed, MPdir)
# install.packages(file.path(MPdir, "pd.Mouse430A2.Mm.ENSG"), repos = NULL, type = "source")

cels <- list.files(file.path(MPdir, "RAW"), pattern = "CEL")

#loading and normalizing data
RawAffy <- read.celfiles(file.path(MPdir, "RAW", cels), pkgname = "pd.Mouse430A2.Mm.ENSG")

#perform RMA
eSet = rma(RawAffy)
colnames(eSet) <- sub(".CEL", "", colnames(eSet))
eSet <- eSet[!grepl("AFFX", rownames(eSet)),]

# Create phenotype data

### Get sample information from cel file names
iLines <- readLines(file.path(MPdir, "GSE22033_series_matrix.txt"))
iLines <- iLines[grepl("Sample_source_name_ch1|Sample_geo_accession", iLines)]
iLines <- gsub("!Sample_source_name_ch1|!Sample_geo_accession|\t", "", iLines)
iLines <- gsub("\"\"", "\"", iLines)
iLines <- gsub("\"", ":", iLines)
iLines <- strsplit(iLines, ":")

pdata <- data.frame(sample = iLines[[1]], info =iLines[[2]], row.names = iLines[[1]])[-1,]
pdata$type <- ifelse(grepl("wt", pdata$info), "wt", "mut")
pdata$treat <- "None"; pdata$treat[grepl("rosi", pdata$info)] <- "ROSI"; pdata$treat[grepl("MRL24", pdata$info)] <- "MRL24"

# Create feature data set
fDat <- fData(eSet)
fDat$ENSEMBL <- sub("_at", "", rownames(fDat))

# Add symbols
toHGNC <- select(mouse430a2.db,
                 key=fDat$ENSEMBL, 
                 columns=c("ENSEMBL", "SYMBOL"),
                 keytype="ENSEMBL")
toHGNC <- toHGNC[!duplicated(toHGNC$ENSEMBL),]
rownames(toHGNC) <- toHGNC$ENSEMBL
toHGNC <- toHGNC[fDat$ENSEMBL,]
if(identical(rownames(toHGNC), fDat$ENSEMBL)) fDat$SYMBOL <- toHGNC$SYMBOL

eSet = to.eSet(exprs(eSet), pdata, fDat)

# Run limma for wt v control
ePhos <- eSet[,eSet$treat == "None"]
ePhos$type <- factor(ePhos$type, levels = c("wt", "mut"))
design <- model.matrix(~type, data = pData(ePhos))
fit <- lmFit(ePhos, design)
fit <- eBayes(fit)
out <- topTable(fit, number = Inf, sort.by = "P")
out <- out[!is.na(out$SYMBOL),]

# Create gene sets
Pup <- c("GSE22033_pS273ko_MEF_Up", "", toupper(out$SYMBOL[out$adj.P.Val < 0.05 & out$logFC > 0]))
Pdown <- c("GSE22033_pS273ko_MEF_Down", "", toupper(out$SYMBOL[out$adj.P.Val < 0.05 & out$logFC < 0]))
Pvec <- c(
  paste(Pup, collapse = "\t"),
  paste(Pdown, collapse = "\t")
)

# Create gmt file
writeLines(Pvec, gmtOut_GSE22033)
