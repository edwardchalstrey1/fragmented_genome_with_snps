plot_it <- function(x, y, se, group, best_scores, dataset_run, filename){
	library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
	df <- data.frame(
      generation = factor(x),
	  metric_score = y, # this is just the averages for each population of fitness and metric score
	  group = factor(group),
	  se = se,
	  best = best_scores # these are the best fitness for each pop, and metric score
	  )
	limits <- aes(ymax = metric_score + se, ymin = metric_score - se)
	p <- ggplot(df, aes(colour = group, y = metric_score, x = generation))
	p <- p + geom_point(aes(y = best)) + geom_crossbar(limits, width = 0.2) + facet_grid(group~., scales = "free_y", space = "fixed")
	ggsave(p, file = paste("~/fragmented_genome_with_snps/arabidopsis_datasets/", dataset_run,"/", filename,".png", sep = ""))
}

st_err <- function(a) sd(a)/sqrt(length(a))