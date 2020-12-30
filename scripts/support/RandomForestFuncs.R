# Importance stats
plotImportance <- function(varImp, n = 20){
  impSub <- data.frame(Importance = varImp[1:n,], Gene = rownames(varImp)[1:n])
  impSub$Gene <- factor(impSub$Gene, levels = rev(impSub$Gene[order(impSub$Importance)]))
  p <- ggplot(impSub, aes(Gene, Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    ggtitle(paste0("Gene Importance (Top ", n, ")")) +
    theme_bw()
  return(p)
}

# Classification stats
getClassStats <- function(res, eSet){
  pDat <- pData(eSet)[, c("Chemical", "RXR_Lig", "PPARg_Lig", "PPARg_Mod")]
  Pred <- predict(res$results, t(exprs(eSet)), "prob")
  Pred$Pred <- c("Yes", "No")[apply(Pred[,1:2], 1, which.max)]
  Pred <- cbind(pDat, Pred)
  Pred <- Pred[order(Pred$No),]
  
  # Pool votes for chemicals
  predChems <- as.data.frame(Pred %>% group_by(Chemical) %>%
                               summarise(Yes = mean(Yes), 
                                         No = mean(No), 
                                         RXR_Lig = unique(RXR_Lig),
                                         PPARg_Lig = unique(PPARg_Lig),
                                         PPARg_Mod = unique(PPARg_Mod)))
  predChems$Pred <- c("Yes", "No")[apply(predChems[, c("Yes", "No")], 1, which.max)]
  predChems <- predChems[order(predChems$No),]
  
  # Create summary table
  ## Performance Stats
  sampStats <- res$sampleResults$SensSpec
  Prev <- mean(res$sampleResults$summary$Label == "Yes")
  sampStats <- c(TruePos = Prev,sampStats, AUC = res$sampleResults$roc$auc)
  
  chemStats <- res$chemicalResults$SensSpec
  PrevChem <- mean(res$chemicalResults$summary$Label == "Yes")
  chemStats <- c(PropYes = PrevChem, chemStats, AUC = res$chemicalResults$roc$auc)
  
  framStats <- rbind(sampStats , chemStats)
  rownames(framStats) <- c("Individual Samples", "Pooled Votes Per Chemical")
  
  # Finally get variable importance
  varImp <- data.frame(MeanDecreaseGini = res$importanceSorted)
  
  return(list(Pred = Pred, predChems = predChems, framStats = framStats, varImp = varImp))
}
