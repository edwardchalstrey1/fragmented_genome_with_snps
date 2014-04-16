reproduction_plot <- function(x_gen, y_fits, se, group, location, filename){
	library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
	df <- data.frame(
		generation = factor(x_gen),
		y_fits = y_fits,
		group = factor(group),
		se = se
		)
	limits <- aes(ymax = y_fits + se, ymin = y_fits - se)
	p <- ggplot(df, aes(colour = group, y = y_fits, x = generation)) + geom_crossbar(limits, width = 0.2)
	ggsave(p, file = paste("~/", location, filename,".png", sep = ""))
}

st_err <- function(a) sd(a)/sqrt(length(a))