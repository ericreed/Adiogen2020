require(ggplot2)
require(pROC)
require(tidyr)
require(dplyr)

# Directories
datDir <- "../data"
resDir <- "../results"
supDir <- "support"

## Subdirectories
outputDir <- file.path(resDir, "CVoutput")
outputTabDir <- file.path(outputDir, "tables")
outputPlotDir <- file.path(outputDir, "plots")

# Files to read in
kCViFiles <- list.files(outputDir)

# Files to write out
cvOutPutFile <- file.path(outputTabDir, "cvPerf.csv")
cvOutStatFile <- file.path(outputTabDir, "cvPerf_stats.csv")
cvOutBoxPlot <- file.path(outputPlotDir, "cvPerfBoxplot.png")
rocBMout <- file.path(outputPlotDir, "begMergeROCcurves.png")

# Read in helper functions
source(file.path(supDir, "summarise_kCV.R"))

# Create file list for each type of run
types <- c("BagMerge", "PreMerge", "Classic")
typeNames <- c("Bag-merge", "Pre-merge", "Classic")
typeList <- lapply(types, function(x) kCViFiles[grepl(x, kCViFiles)])
names(typeList) <- types

# Get results
resList <- do.call(rbind, lapply(types, function(x) {
    typeFiles <- typeList[[x]]
    typeFram <- do.call(rbind, lapply(typeFiles, function(y, x) {
        res <- summarise_kCVi(readRDS(file.path(outputDir, y)), collapse = FALSE)

        # Create data frame of performance
        perf <- res$perfFram[, c("auc", "balacc", "f1", "prec", "sens", "spec")]
        perf$iter <- as.numeric(gsub(".rds", "", strsplit(y, "_")[[1]][length(strsplit(y, "_")[[1]])]))
        perf$type <- typeNames[types == x]

        # Create data frame of predictions
        pred <- res$predFram[, c("votes", "class")]
        pred$iter <- as.numeric(gsub(".rds", "", strsplit(y, "_")[[1]][length(strsplit(y, "_")[[1]])]))
        pred$type <- typeNames[types == x]

        return(list(perf, pred))
    }, x))

    if(x == "Classic") {
        pmFram <- do.call(rbind, lapply(typeFiles, function(y, x) {
            res <- summarise_kCVi(readRDS(file.path(outputDir, y)), collapse = TRUE)

            # Create data frame of performance
            perf <- res$perfFram[, c("auc", "balacc", "f1", "prec", "sens", "spec")]
            perf$iter <- as.numeric(gsub(".rds", "", strsplit(y, "_")[[1]][length(strsplit(y, "_")[[1]])]))
            perf$type <- "Pooled"

            # Create data frame of predictions
            pred <- res$predFram[, c("votes", "class")]
            pred$iter <- as.numeric(gsub(".rds", "", strsplit(y, "_")[[1]][length(strsplit(y, "_")[[1]])]))
            pred$type <- "Pooled"

            return(list(perf, pred))
        }, x))
        typeFram <- rbind(pmFram, typeFram)
    }
    kCVlist <- list(do.call(rbind, typeFram[, 1]),
                    do.call(rbind, typeFram[, 2]))

    return(kCVlist)
}))

# Get performance and rename columns
resFram <- do.call(rbind, resList[,1])
resFram <- resFram[, c("type", "iter", "auc", "balacc", "f1", "prec", "sens", "spec")]
colnames(resFram) <- c("type", "iter", "AUC", "Bal. Acc.", "F1", "Precision", "Sensitivity", "Specificity")

# Save table
write.csv(resFram, cvOutPutFile, row.names = FALSE)

# Reformat data for boxplots
resGath <- resFram %>%
    gather("metric", "value", c("AUC", "Bal. Acc.", "F1", "Precision", "Sensitivity", "Specificity"))
resGath$type<- factor(resGath$type, levels = c("Bag-merge", "Pre-merge", "Pooled", "Classic"))

# Create boxplots
png(cvOutBoxPlot, height = 350, width = 900)
ggplot(resGath, aes(y = value, x = type)) +
    geom_boxplot() +
    facet_wrap(~metric, nrow = 1) +
    scale_y_continuous(name = "Performance Estimate") +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold", size = 18),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(face = "bold", size = 15),
        panel.spacing = unit(2, "lines")
    )
dev.off()

# Get mean and sd estimates
resSum <- resGath %>%
    group_by(type, metric) %>%
    summarise(
        mean = mean(value),
        sd = sd(value)
    )
write.csv(resSum, cvOutStatFile, row.names = FALSE)

# Plot ROC curves for bag merge

## Get predictions
predBM <- do.call(rbind, resList[,2])
predBM <- predBM[predBM$type == "Bag-merge",]

## Add sensitivity and specificity
ssBM <- do.call(rbind, lapply(unique(predBM$iter), function(iter) {

    # Subset for iteration
    iterFram <- predBM[predBM$iter == iter,]
    iterFram <- iterFram[order(iterFram$votes),]

    # Get sens and spec at each threshold
    statFram <- do.call(rbind, lapply(2:nrow(iterFram), function(x, iterFram){
        sens <- sum(iterFram$class[x:nrow(iterFram)] == "Yes") / sum(iterFram$class == "Yes")
        spec <- sum(iterFram$class[1:(x-1)] == "No") / sum(iterFram$class == "No")
        data.frame(sens = sens, spec = spec)
    }, iterFram))
    statFram <- statFram[!duplicated(statFram$sens),]
    statFram <- statFram[!duplicated(statFram$spec),]
    statFram <- rbind(data.frame(sens = 1, spec = 0),
                      statFram,
                      data.frame(sens = 0, spec = 1))
    statFram$iter <- iter
    return(statFram)
}))

## Generate step plot
png(rocBMout)
ggplot(ssBM, aes(x = spec, y = sens, group = as.factor(iter))) +
    geom_step(direction = "vh", alpha = 0.4) +
    scale_y_continuous(name = "Sensitivity") +
    scale_x_continuous(name = "Specificity") +
    geom_abline(intercept = 1, slope = -1, color = "red") +
    theme_bw() +
    theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
    )
dev.off()
