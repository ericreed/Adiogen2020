plot_cv_summary <- function(cv_summary){

  sumFram <- as.data.frame( do.call(rbind, lapply(names(cv_summary), function(summName){

    # Get summarised results
    summ <- cv_summary[[summName]]

    # Add rep id to results
    for(i in 1:length(summ)) summ[[i]]$sub_tree_performance$repID <- as.character(i)

    # Combine reps
    cvSum <- do.call(rbind, lapply(summ, function(x) x$sub_tree_performance))

    # Gather statistics
    statVec <- c("auc", "balacc", "f1", "prec", "sens", "spec", "acc")
    cvGath <- cvSum %>% gather(key = "metric", value = "stat", statVec)
    cvGath$metric <- factor(cvGath$metric, levels = statVec)

    # Add method name
    cvGath$ID <- summName

    return(cvGath)
  })))
  sumFram$ID <- factor(sumFram$ID, levels = names(cv_summary))

  # Set breaks for x-axis
  maxTrees <- max(sumFram$fTrees)
  minTrees <- min(sumFram$fTrees)
  breaksX <- seq(0, maxTrees, by = 25)

  # Add values at max trees
  sumFram$maxStat <- sumFram$stat;   sumFram$maxStat[sumFram$fTrees < maxTrees] <- NA

  # Force two decimal places
  scaleFUN <- function(x) sprintf("%.2f", x)

  # Plot results in line plot
  pLines <- suppressWarnings(

    ggplot(sumFram, aes(x = fTrees, y = stat)) +
      geom_smooth(se = F, size = 1.5) +
      geom_line(aes(group = repID), alpha = 0.3) +
      geom_boxplot(aes(x = maxTrees, y = maxStat), width = 0.5) +
      scale_x_continuous(name = "Number of Trees", breaks = breaksX) +
      scale_y_continuous(name = "", labels = scaleFUN) +
      facet_grid(metric~ID, scales = "free_y")

  )

  # Simple boxplots
  pBox <- suppressWarnings(

    ggplot(sumFram, aes(x = ID, y = maxStat)) +
      geom_boxplot(width = 0.5) +
      scale_y_continuous(name = "", labels = scaleFUN) +
      facet_grid(~metric)

  )

  return(list(lineplot = pLines, boxplot = pBox))
}
