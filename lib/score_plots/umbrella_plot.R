av_sd <- function(df, metric){
  names(df)[names(df)==metric] <- "metric_scores"
  df <- transform(df, Average = ave(metric_scores, gen, replicates, param_types))
  sd = aggregate(df$metric_scores, by=list(gen=df$gen, replicates=df$replicates, param_types=df$param_types), sd)
  df <- merge(df,sd)
  return(df)
}


uplot <- function(df, title){
	library(ggplot2)
	p <- ggplot(df, aes(colour = replicates, y = metric_scores, x = gen)) +
	    xlab("Generations") +
	    ylab(metric) +
	    ggtitle(title) +
	    scale_y_continuous(limits=c(0, 1)) +
		scale_x_continuous() +
		geom_line(aes(y = Average), size=1) +
		geom_line(aes(y = Average-x), size=0.5, linetype="dotted") +
		geom_line(aes(y = Average+x), size=0.5, linetype="dotted") +
	    facet_grid(param_types~., scales = "free_y", space = "fixed") +
	    theme_bw()
	return(p)
}