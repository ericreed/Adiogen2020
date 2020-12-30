require(openxlsx)
require(ggplot2)
require(Biobase)
require(dplyr)
source("parameters.R")

# Set script tag
sTag <- "779b"

# Begin

# Read in file to update chemical names
update <- read.csv(updateFile, stringsAsFactors = F)
update$Abbreviation <- toupper(update$Abbreviation)
update$Abbreviation <- sub("-", "_", update$Abbreviation)

# Subset for just updates
update <- update[, c("Abbreviation", "Update")]
colnames(update) <- c("Chemical", "New")

# Read in expression set
eSet <- readRDS(eSetCombatFile011a)

# Rename vehicle
eSet$PPARg_Mod[grepl("VH", eSet$Chemical)] <- "Vehicle"
eSet$Chemical[grepl("VH", eSet$Chemical)] <- "Vehicle"

# Rename naive
eSet$PPARg_Mod[eSet$Chemical == "NAIVE"] <- "Naive"


# Read in random forest results
gImp <- read.xlsx(excelFileOut779a, sheet = "Gene_Importance", rowNames = T)

# Get top two genes
topGenes <- row.names(gImp)[seq(2)]

# Function to plot genes
plotGenes <- function(gene) {
  
  # Get data.frame of expression and chemical information
  df <- data.frame(t(exprs(eSet[gene,])), pData(eSet)[, c("Chemical", "PPARg_Mod")], stringsAsFactors = FALSE)
  df$`Sample ID` <- rownames(df)
  colnames(df)[1] <- "exp"
  
  # Merge with update
  df <- merge(df, update, all = TRUE)
  df$New[is.na(df$New)] <- df$Chemical[is.na(df$New)]
  
  # Combine Veihcle and Naive for coloring
  df$PPARg_Mod_VN <- df$PPARg_Mod
  df$PPARg_Mod_VN[df$PPARg_Mod_VN %in% c("Vehicle", "Naive")] <- "Vehicle/Naive"
  df$PPARg_Mod_VN <- factor(df$PPARg_Mod_VN, levels = c("Yes", "Suspected", "No", "Vehicle/Naive"))
  
  # Get means for each chemical
  dfMeans <- df %>%
    group_by(New, PPARg_Mod_VN) %>%
    summarise(Mexp = mean(exp))
  
  # Order chemicals by mean
  dfMeans <- dfMeans[order(dfMeans$Mexp, decreasing = TRUE),]
  dfMeans$New <- factor(dfMeans$New, levels = dfMeans$New)
  df$New <- factor(df$New, levels = levels(dfMeans$New))
  
  # Plot
  p <- ggplot(data = df, aes(x = New, y = exp, group = New, color = PPARg_Mod_VN, shape = PPARg_Mod_VN)) +
    geom_hline(yintercept = unlist(dfMeans[dfMeans$New == "Vehicle", "Mexp"]), color = "grey70", size = 0.8) +
    geom_point(data = dfMeans, aes(x = New, y = Mexp), shape = 3, size = 3, color = "black") +
    geom_line(color = "grey60", size = 1) +
    geom_point(size = 3) +
    scale_color_manual(name = "PPARγ Activity Modifer",
                       values = c(
                         "Vehicle/Naive" = "grey50",
                         "Yes" = "red",
                         "No" = "blue",
                         "Suspected" = "chartreuse4"
                       )) +
    scale_shape_manual(name = "PPARγ Activity Modifer",
                       values = c(
                         "Vehicle/Naive" = 18,
                         "Yes" = 16,
                         "No" = 17,
                         "Suspected" = 15
                       )) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
      axis.text.y = element_text(size = 13),
      legend.position="bottom",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20)
    )
  
  png(paste0(GenePlotPrefix779b, gene, "_", sTag, ".png"), width = 1200, height = 400)
  print(p)
  dev.off()
  
  # Return data frame
  df$PPARg_Mod[df$PPARg_Mod == "Vehicle"] <- "No"
  df <- df[, c("Sample ID", "New", "PPARg_Mod", "exp")]
  colnames(df) <- c("Sample ID", "Chemical", "class", gene)
  return(df)
}


# Get data frame of sample expression levels
dfFram <- Reduce(merge, lapply(topGenes, plotGenes))

# Order by first gene
dfFram <- dfFram[order(dfFram[, topGenes[1]], decreasing = TRUE),]
dfTopMeans <- dfFram %>%
  group_by(Chemical) %>%
  summarise(mExp = mean(eval(parse(text = topGenes[1]))))
dfTopMeans <- dfTopMeans[order(dfTopMeans$mExp, decreasing = TRUE),]
dfFram <- dfFram[order(match(dfFram$Chemical, dfTopMeans$Chemical)),]

# Add to workbook
wb <- loadWorkbook(excelFileOut779a)
addWorksheet(wb, "Expression Top Genes")
writeData(wb, "Expression Top Genes", dfFram, startRow = 8)
saveWorkbook(wb, excelFileOut779b, overwrite = TRUE)
