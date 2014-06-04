uplot <- function(generations, metric_scores, run){
  library(ggplot2)
  df <- data.frame(
    x = factor(generations),
    y = metric_scores,
    group = factor(run)
    )
  ggplot(df, aes(colour = group, y = y, x = x)) +
    geom_point(aes(y = y)) +
    geom_smooth(aes(y = y, group=group), method="loess") # choose most appropriate smoothing method 
}

mets <- c(0.55, 0.60, 0.70, 0.65, 0.70, 0.80, 0.30, 0.40, 0.50, 0.40, 0.50, 0.60)
gens <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2)
groups <- c("run1", "run1", "run1", "run1", "run1", "run1", "run2", "run2", "run2", "run2", "run2", "run2")

uplot(gens,mets,groups)