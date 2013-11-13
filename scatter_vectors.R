scatta <- function(y, method_name, method_id, dataset){
x <- (1:length(y))
title <- paste("Graph to show SNP densities for fragments ordered by: ", method_name, sep="")
png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/", dataset, "/figures/", method_id, ".png", sep=""))
graph <- (plot(x, y, main=title, xlab="Fragments", ylab="SNP Density (SNPs/Kb)"))
lines(lowess(x, y), col="red", lwd=4)
dev.off()
return(graph)
}
