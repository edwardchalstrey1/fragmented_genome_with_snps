uplot <- function(generations, metric_scores, runs){
	library(ggplot2)
	df <- data.frame(
		x = factor(generations),
		y = metric_scores,
		group = factor(runs)
    )
	p <- ggplot(df, aes(colour = group, y = y, x = x))
	p <- p + geom_point(aes(y = y))
	return(p)
}

mets <- c(0.55, 0.60, 0.70, 0.65, 0.70, 0.80, 0.30, 0.40, 0.50, 0.40, 0.50, 0.60)
gens <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2)
groups <- c("run1", "run1", "run1", "run1", "run1", "run1", "run2", "run2", "run2", "run2", "run2", "run2")

uplot(gens,mets,groups)