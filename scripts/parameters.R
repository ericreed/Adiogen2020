# Main directories
baseDir <- ".."
supDir <- "support"
datDir <- file.path(baseDir, "data")
resDir <- file.path(baseDir, "results")
plotDir <- file.path(resDir, "plots")
tabDir <- file.path(resDir, "tables")
clasDir <- file.path(resDir, "Classification")
taxDir <- file.path(resDir, "Taxonomy")
packDir <- "../HierForest/"

## Subdirectories
esetDir <- file.path(datDir, "eSets")
gsigDir <- file.path(datDir, "GenSigs")
gsetDir <- file.path(datDir, "GeneSets")
datMatDir <- file.path(datDir, "dataMatrix")
rdatDir <- file.path(resDir, "RDATA")
ppDir <- file.path(plotDir, "PreProc")
boxDir037a <- file.path(plotDir, paste0("DGE_ANOVA_boxplots_037a"))
boxDir98 <- file.path(plotDir, paste0("ssGSEA_boxplots_098a"))
boxDir111a <- file.path(plotDir, paste0("ssGSEA_boxplots_111a"))
boxDir508a <- file.path(plotDir, paste0("ssGSEA_boxplots_508a"))
nrPlotDir020a <- file.path(plotDir, paste0("RegressNileRed_MixedModel_020a"))
graphDir <- file.path(plotDir, "graphs")
classTab <- file.path(clasDir, "tables")
classPlot <- file.path(clasDir, "plots")
taxTab <- file.path(taxDir, "tables")
taxPlot <- file.path(taxDir, "plots")
curGS <- file.path(datDir, "Curated Genesets")
chemDir <- file.path(datDir, "ChemInfo")
sampDir <- file.path(datDir, "SampleInfo")
nrDir <- file.path(datDir, "NileRed")
MPdir <- file.path(datDir, "GSE22033")
HMAdir <- file.path(datDir, "Human_28257690")

# R function files
source(file.path(supDir, "helper.R"))
source(file.path(supDir, "CreateGeneSignature.R"))
source(file.path(supDir, "RandomForestFuncs.R"))
source(file.path(supDir, "TaxonomyFunctions.R"))

# Main files
goFile <- file.path(gsetDir, "c5.bp.v6.2.symbols.gmt")
SampFile <- file.path(sampDir, "SampleInfo.csv")
updateFile <- file.path(chemDir, "ChemInfo.csv")

# 000
eSetFile000 <- file.path(esetDir, paste0("Obesogen_3DGE_", "plate1", "_", "000", ".rds"))
eSetRemFile000 <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_", "plate1", "_", "000", ".rds"))
AboxplotFile000 <- file.path(ppDir, paste0("Aboxplot_", "plate1", "_", "000", ".png"))
BboxplotFile000 <- file.path(ppDir, paste0("Bboxplot_", "plate1", "_", "000", ".png"))
AbarplotFile000 <- file.path(ppDir, paste0("Abarplot_", "plate1", "_", "000", ".png"))
BbarplotFile000 <- file.path(ppDir, paste0("Bbarplot_", "plate1", "_", "000", ".png"))
medVmadplotFile000 <- file.path(ppDir, paste0("medVmadplot_", "plate1", "_", "000", ".png"))
clustMADplotFile000 <- file.path(ppDir, paste0("clustMADplot_", "plate1", "_", "000", ".png"))
clustExpplotFile000 <- file.path(ppDir, paste0("clustExpplot_", "plate1", "_", "000", ".png"))

# 001
eSetFile001 <- file.path(esetDir, paste0("Obesogen_3DGE_", "plate2", "_", "001", ".rds"))
eSetRemFile001 <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_", "plate2", "_", "001", ".rds"))
AboxplotFile001 <- file.path(ppDir, paste0("Aboxplot_", "plate2", "_", "001", ".png"))
BboxplotFile001 <- file.path(ppDir, paste0("Bboxplot_", "plate2", "_", "001", ".png"))
AbarplotFile001 <- file.path(ppDir, paste0("Abarplot_", "plate2", "_", "001", ".png"))
BbarplotFile001 <- file.path(ppDir, paste0("Bbarplot_", "plate2", "_", "001", ".png"))
medVmadplotFile001 <- file.path(ppDir, paste0("medVmadplot_", "plate2", "_", "001", ".png"))
clustMADplotFile001 <- file.path(ppDir, paste0("clustMADplot_", "plate2", "_", "001", ".png"))
clustExpplotFile001 <- file.path(ppDir, paste0("clustExpplot_", "plate2", "_", "001", ".png"))

# 002
eSetFile002 <- file.path(esetDir, paste0("Obesogen_3DGE_", "plate3", "_", "002", ".rds"))
eSetRemFile002 <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_", "plate3", "_", "002", ".rds"))
AboxplotFile002 <- file.path(ppDir, paste0("Aboxplot_", "plate3", "_", "002", ".png"))
BboxplotFile002 <- file.path(ppDir, paste0("Bboxplot_", "plate3", "_", "002", ".png"))
AbarplotFile002 <- file.path(ppDir, paste0("Abarplot_", "plate3", "_", "002", ".png"))
BbarplotFile002 <- file.path(ppDir, paste0("Bbarplot_", "plate3", "_", "002", ".png"))
medVmadplotFile002 <- file.path(ppDir, paste0("medVmadplot_", "plate3", "_", "002", ".png"))
clustMADplotFile002 <- file.path(ppDir, paste0("clustMADplot_", "plate3", "_", "002", ".png"))
clustExpplotFile002 <- file.path(ppDir, paste0("clustExpplot_", "plate3", "_", "002", ".png"))

# 003
eSetFile003 <- file.path(esetDir, paste0("Obesogen_3DGE_", "plate4", "_", "003", ".rds"))
eSetRemFile003 <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_", "plate4", "_", "003", ".rds"))
AboxplotFile003 <- file.path(ppDir, paste0("Aboxplot_", "plate4", "_", "003", ".png"))
BboxplotFile003 <- file.path(ppDir, paste0("Bboxplot_", "plate4", "_", "003", ".png"))
AbarplotFile003 <- file.path(ppDir, paste0("Abarplot_", "plate4", "_", "003", ".png"))
BbarplotFile003 <- file.path(ppDir, paste0("Bbarplot_", "plate4", "_", "003", ".png"))
medVmadplotFile003 <- file.path(ppDir, paste0("medVmadplot_", "plate4", "_", "003", ".png"))
clustMADplotFile003 <- file.path(ppDir, paste0("clustMADplot_", "plate4", "_", "003", ".png"))
clustExpplotFile003 <- file.path(ppDir, paste0("clustExpplot_", "plate4", "_", "003", ".png"))

# 004
eSetFile004 <- file.path(esetDir, paste0("Obesogen_3DGE_", "plate5", "_", "004", ".rds"))
eSetRemFile004 <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_", "plate5", "_", "004", ".rds"))
AboxplotFile004 <- file.path(ppDir, paste0("Aboxplot_", "plate5", "_", "004", ".png"))
BboxplotFile004 <- file.path(ppDir, paste0("Bboxplot_", "plate5", "_", "004", ".png"))
AbarplotFile004 <- file.path(ppDir, paste0("Abarplot_", "plate5", "_", "004", ".png"))
BbarplotFile004 <- file.path(ppDir, paste0("Bbarplot_", "plate5", "_", "004", ".png"))
medVmadplotFile004 <- file.path(ppDir, paste0("medVmadplot_", "plate5", "_", "004", ".png"))
clustMADplotFile004 <- file.path(ppDir, paste0("clustMADplot_", "plate5", "_", "004", ".png"))
clustExpplotFile004 <- file.path(ppDir, paste0("clustExpplot_", "plate5", "_", "004", ".png"))

# 005
eSetFile005 <- file.path(esetDir, paste0("Obesogen_3DGE_", "plate6", "_", "005", ".rds"))
eSetRemFile005 <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_", "plate6", "_", "005", ".rds"))
AboxplotFile005 <- file.path(ppDir, paste0("Aboxplot_", "plate6", "_", "005", ".png"))
BboxplotFile005 <- file.path(ppDir, paste0("Bboxplot_", "plate6", "_", "005", ".png"))
AbarplotFile005 <- file.path(ppDir, paste0("Abarplot_", "plate6", "_", "005", ".png"))
BbarplotFile005 <- file.path(ppDir, paste0("Bbarplot_", "plate6", "_", "005", ".png"))
medVmadplotFile005 <- file.path(ppDir, paste0("medVmadplot_", "plate6", "_", "005", ".png"))
clustMADplotFile005 <- file.path(ppDir, paste0("clustMADplot_", "plate6", "_", "005", ".png"))
clustExpplotFile005 <- file.path(ppDir, paste0("clustExpplot_", "plate6", "_", "005", ".png"))

# 010a
eSetRawFile010a <- file.path(esetDir, paste0("Obesogen_3DGE_LERem_raw_", "010a", ".rds"))
HM_CPM_MADFile010a <- file.path(ppDir, paste0("Heatmap_allPlates_LERem_logCPM_bothExperiments_3Kmad_", "010a", ".png"))
HM_CPM_MEDFile010a <- file.path(ppDir, paste0("Heatmap_allPlates_LERem_logCPM_bothExperiments_3KmedExp_", "010a", ".png"))
eSetAFile010a <- file.path(esetDir, paste0("Obesogen_3DGE_ExpA_raw_LERem_", "010a", ".rds"))
eSetBFile010a <- file.path(esetDir, paste0("Obesogen_3DGE_ExpB_raw_LERem_", "010a", ".rds"))

LibrarySizeAllChems <- file.path(ppDir, paste0("BoxPlots_LibrarySizes_AllChems_ByPlate_ExpA", "010a", ".png"))
LibrarySizeAllChems_NT <- file.path(ppDir, paste0("BoxPlots_LibrarySizes_AllChems_NotTransformed_ByPlate_ExpA", "010a", ".png"))
LibrarySizeByChem <- file.path(ppDir, paste0("BoxPlots_LibrarySizes_AllChems_ByChem_ExpA", "010a", ".png"))
LibrarySizeVehicles <- file.path(ppDir, paste0("BoxPlots_LibrarySizes_VehiclesOnly_ByPlate_ExpA", "010a", ".png"))
LibrarySizePlate5 <- file.path(ppDir, paste0("BoxPlots_LibrarySizes_Plate5_ByChem_ExpA", "010a", ".png"))
LibrarySizePlate6 <- file.path(ppDir, paste0("BoxPlots_LibrarySizes_Plate6_ByChem_ExpA", "010a", ".png"))

# 011a
MADvExpFile011a <- file.path(ppDir, paste0("ExpA_MADvExp_LERem_", "011a", ".png"))
combatFiles011a <- file.path(ppDir, paste0("ExpA_Combat_LERem_"))
Cluster_Prefilter1 <- file.path(ppDir, paste0("Sample_Dendro_combat_1_", "011a", ".png"))
Cluster_RemoveSamps1 <- file.path(ppDir, paste0("Sample_Dendro_combat_1_remSamps", "011a", ".csv"))
MADvExpFile_2_011a <- file.path(ppDir, paste0("ExpA_MADvExp_LERem_PostFiler1_", "011a", ".png"))
combatFiles_2_011a <- file.path(ppDir, paste0("ExpA_Combat_LERem_PostFiler1_"))
Cluster_Prefilter2 <- file.path(ppDir, paste0("Sample_Dendro_combat_2_", "011a", ".png"))
Cluster_RemoveSamps2 <- file.path(ppDir, paste0("Sample_Dendro_combat_2_remSamps", "011a", ".csv"))
MADvExpFile_3_011a <- file.path(ppDir, paste0("ExpA_MADvExp_LERem_PostFiler2_", "011a", ".png"))
combatFiles_3_011a <- file.path(ppDir, paste0("ExpA_Combat_LERem_PostFiler2_"))
Cluster_Prefilter3 <- file.path(ppDir, paste0("Sample_Dendro_combat_3_", "011a", ".png"))
Cluster_RemoveSamps3 <- file.path(ppDir, paste0("Sample_Dendro_combat_3_remSamps", "011a", ".csv"))
MADvExpFile_4_011a <- file.path(ppDir, paste0("ExpA_MADvExp_LERem_PostFiler3_", "011a", ".png"))
combatFiles_4_011a <- file.path(ppDir, paste0("ExpA_Combat_LERem_PostFiler3_"))
Cluster_Prefilter4 <- file.path(ppDir, paste0("Sample_Dendro_combat_4_", "011a", ".png"))
eSetCombatFile011a <- file.path(esetDir, paste0("Obesogen_3DGE_ExpA_AllSamples_Combat_LERem_", "011a", ".rds"))

# Differential Analysis

# 032a
LIMMAoutFile032a <- file.path(tabDir, paste0("UnivariateDiffAnalysis_ExpA_LowExpRem_Combat_LERem_", "032a", ".xlsx"))
LIMMAoutFileSig032a <- file.path(tabDir, paste0("UnivariateDiffAnalysis_FDR1_ExpA_LowExpRem_Combat_LERem_", "032a", ".xlsx"))
dataMatrixFile032a <- file.path(datMatDir, paste0("UnivariateDiffAnalysis_testStatistics_ExpA_LowExpRem_Combat_LERem_", "032a", ".rds"))

# Random Forest Classification

## 779a
ForestOut779a <- file.path(clasDir, paste0("RandomForest_ExpA_LowExpRem_BagMerge", "779a", ".rds"))
excelFileOut779a <- file.path(classTab, paste0("RandomForestClassification_ExpA_LowExpRem_BagMerge", "779a", ".xlsx"))

## 779b
GenePlotPrefix779b <- file.path(classPlot, "RandomForest_TopGene_")
excelFileOut779b <- file.path(classTab, paste0("RandomForestClassification_ExpA_LowExpRem_BagMerge_TopGenes_", "779b", ".xlsx"))


# K2Taxonomer

# 893a
SplitSummaryTableOut <- file.path(taxTab, paste0("Taxonomy_Results_", "893a", ".xlsx"))
dendroPlotOut <- file.path(taxPlot, paste0("TaxonomyDendrogram_", "893a", ".png"))
HclustPlotOut <- file.path(taxPlot, paste0("HclustDendrogram_", "893a", ".png"))
rData893a <- file.path(taxDir, paste0("TaxRes_", "893a", ".RData"))

# Integrative data analysis

# s273 KO data
gmtOut_GSE22033 <- file.path(gsetDir, "GSE22033_pS273ko.gmt")

# Human Measurements data
varInfoFile <- file.path(HMAdir, "VariableInfo.csv")
Human_key_file <- file.path(HMAdir, "key.png")
eSet_Human_124 <- file.path(esetDir, "Human_28257690_eSet_124.rds")
eSet_Human_ssGSEA_tax_124 <- file.path(esetDir, "Human_28257690_eSet_ssGSEA_tax_124.rds")
Human_ssGSEA_tax_res_124 <- file.path(tabDir, "Human_28257690_eSet_ssGSEA_tax_124.xlsx")
Human_ssGSEA_Heatmap <- file.path(plotDir, "dendros_124", "Human_ssGSEA_tax_res_124_heatmap.png")
Human_ssGSEA_Lines <- file.path(plotDir, "dendros_124", "Human_ssGSEA_tax_res_124_lines.png")
Human_ssGSEA_Lines2 <- file.path(plotDir, "dendros_124", "Human_ssGSEA_tax_res_124_lines2.png")
Human_ssGSEA_tax_res_prefix <- file.path(plotDir, "dendros_124", "Human_ssGSEA_tax_res_124_")
