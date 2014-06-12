uplot <- function(generations, metric_scores, runs, param_types, metric, title){
	library(ggplot2)
	df <- data.frame(
		gen = generations,
		metric_scores = metric_scores,
		replicates = runs,
    param_types = param_types
    )
	p <- ggplot(df, aes(colour = replicates, y = metric_scores, x = gen)) +
    xlab("Generations") +
    ylab(metric) +
    ggtitle(title) +
	  # geom_point(aes(y = metric_scores,x = factor(gen))) +
	  geom_boxplot(aes(y = metric_scores, x = factor(gen), factor(replicates))) +
    scale_y_continuous(limits=c(0, 1)) +
    scale_x_discrete() +
    facet_grid(param_types~., scales = "free_y", space = "fixed")
  return(p)
}

mets <- c(0.55, 0.60, 0.70, 0.65, 0.70, 0.80, 0.30, 0.40, 0.50, 0.40, 0.50, 0.60)
mets <- c(mets, mets)
gens <- c(1, 1, 1, 10, 10, 10, 1, 1, 1, 10, 10, 10)
gens <- c(gens, gens)
runs <- c("run1", "run1", "run1", "run1", "run1", "run1", "run2", "run2", "run2", "run2", "run2", "run2")
runs <- c(runs, runs)
pt <- c('p1','p1','p1','p1','p1','p1','p1','p1','p1','p1','p1','p1','p2','p2','p2','p2','p2','p2','p2','p2','p2','p2','p2','p2')

uplot(gens,mets,runs,pt,'Fitness', 'Multiple replicates of genetic algorithm to re-order contigs...')