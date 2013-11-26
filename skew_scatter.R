scatta <- function(gradients, nu_ids, filename, grad_string, min_snps_string){
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/skew_scatter_", filename, ".png", sep=""))
	plot(nu_ids, gradients, main=paste("The ", grad_string, " of SNPs 
			for fragments with ", min_snps_string, " or more SNPs from dataset 5", sep=""), xlab="Fragments", ylab="Gradient")
	lines(lowess(nu_ids, gradients), col="red", lwd=4)
	dev.off()
}

example <- function(frag_num_string, x, y, min_snps_string){
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/example_gradient_f", frag_num_string, "_mins", min_snps_string, ".png", sep=""))
	plot(x, y, main=paste("Example of gradient determination for fragment ", frag_num_string, " of dataset 5
		(minimum number of SNPs per fragment: ", min_snps_string, ")", sep=""), xlab="Fragment (bases)", ylab="SNP density")
	abline(coef=coef(lm(y~x)), col="blue", lwd=4)
	dev.off()
}