av_sd <- function(df, metric){
  names(df)[names(df)==metric] <- "metric_scores"
  df <- transform(df, Average = ave(metric_scores, gen, replicates, param_types))
  sd = aggregate(df$metric_scores, by=list(gen=df$gen, replicates=df$replicates, param_types=df$param_types), sd)
  df <- merge(df,sd)
  return(df)
}

uplot <- function(df, title, y_axis){
	df$param_types = factor(df$param_types, levels=c('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12'))
	df$replicates<- as.numeric(sub("p_run","",df[,2]))
	library(ggplot2)
	p <- ggplot(df, aes(colour = factor(replicates), y = metric_scores, x = gen)) +
	    xlab("Generations of genetic algorithm - each generation is a new population of contig permutations") +
	    ylab(y_axis) +
	    ggtitle(title) +
	    scale_y_continuous(limits=c(0, 1)) +
		scale_x_continuous() +
		geom_line(aes(y = 1.0-(Average)), size=0.5) +
		geom_line(aes(y = 1.0-(Average-x)), size=0.1, linetype="solid") +
		geom_line(aes(y = 1.0-(Average+x)), size=0.1, linetype="solid") +
    	# geom_ribbon(aes(y = Average, ymin = (Average-x), ymax = (Average+x), fill = replicates, alpha = 0.0)) +
    	facet_wrap(~param_types, ncol=3) +
	    theme_bw() +
	    theme(title = element_text(size = rel(0.5))) +
		guides(col = guide_legend(keywidth = 0.25, keyheight = 0.25, ncol = 4, byrow = TRUE, title.theme = element_text(size=8, angle = 0),
                              title = "Replicates"))
	return(p)
}

av_sd2 <- function(df, metric){
  names(df)[names(df)==metric] <- "metric_scores"
  df <- transform(df, Average = ave(metric_scores, population))
  sd = aggregate(df$metric_scores, by=list(population=df$population), sd)
  df <- merge(df,sd)
  return(df)
}

metric_test_plot <- function(df, title, y_axis, metric){
	df <- av_sd2(df, metric)
	library(ggplot2)
	p <- ggplot(df, aes(y = metric_scores, x = population)) +
	    xlab('Number of adjacent swaps') + # TODO improve axis names and title
	    ylab(y_axis) +
	    ggtitle(title) +
	    scale_y_continuous(limits=c(0, 1.2)) +
		scale_x_continuous() +
	  geom_line(aes(y = Average), size=0.5) +
	  geom_ribbon(aes(y = Average, ymin = (Average-x), ymax = (Average+x), fill="red", alpha=0.1)) +
	    theme_bw()
  return(p)
}