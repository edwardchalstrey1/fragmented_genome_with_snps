plot_it <- function(x, y, se, group, dataset_run){
	library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
	df <- data.frame(
      generation = factor(x),
	  metric_score = y,
	  group = factor(group),
	  se = se
	  )
	limits <- aes(ymax = metric_score + se, ymin=metric_score - se)
	p <- ggplot(df, aes(colour=group, y=metric_score, x=generation))
	p <- p + geom_point() + geom_errorbar(limits, width=0.2)
	ggsave(p, file = paste("~/fragmented_genome_with_snps/arabidopsis_datasets/", dataset_run,"/met_scores.png", sep=""))
}

st_err <- function(a) sd(a)/sqrt(length(a))