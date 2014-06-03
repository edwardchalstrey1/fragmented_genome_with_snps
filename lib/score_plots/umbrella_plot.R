uplot <- function(generations, metric_scores, group){
  library(ggplot2)
  df <- data.frame(
    x = factor(generations),
    y = metric_scores,
    group = factor(group)
    )
  ggplot(df, aes(colour = group, y = y, x = x)) +
    geom_point(aes(y = y)) +
    geom_smooth(aes(y = y), method="loess") # choose most appropriate smoothing method 
}