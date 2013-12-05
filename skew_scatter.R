scatta <- function(gradients, nu_ids, filename, grad_string, min_snps_string){
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/skew_scatter_", filename, ".png", sep=""))
	plot(nu_ids, gradients, main=paste("The ", grad_string, " of SNPs 
			for fragments with ", min_snps_string, " or more SNPs from dataset 5", sep=""), xlab="Fragments", ylab="Gradient")
	lines(lowess(nu_ids, gradients), col="red", lwd=4)
	dev.off()
}

example <- function(frag_num_string, x, y, min_snps_string, d, m){
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/example_gradient_f", frag_num_string, "_mins", min_snps_string, "_", d, "_m", m, ".png", sep=""))
	plot(x, y, main=paste("Example of gradient determination for fragment ", frag_num_string, " of dataset 5
		(minimum number of SNPs per fragment: ", min_snps_string, ")", sep=""), xlab="Fragment (bases)", ylab="SNP density")
	abline(coef=coef(lm(y~x)), col="blue", lwd=4)
	dev.off()
}

how_scatta <- function(frag_num_string, snp_pos, y, length){
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/how_scatta_f", frag_num_string, ".png", sep=""))
	plot(snp_pos, y, main=paste("Example of SNP distribution for fragment ", frag_num_string, " of dataset5
		(fragment length = ", length, ")", sep=""), xlab="SNP positions", ylab="SNP = 1")
	dev.off()
}
gg_hist <- function(frag_num_string, snp_pos, length){
	#png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/gghist_f", frag_num_string, ".png", sep=""))
	names(snp_pos) <- c("pos")
	g <- ggplot(snp_pos, aes(pos)) + geom_histogram() + ggtitle(quote(paste("Fragment ", frag_num_string, " of dataset5
		(fragment length = ", length, ")" sep=""))) + scale_y_continuous() + scale_x_continuous(name=quote("SNP positions"))
	ggsave(g, file=paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset5/figures/gghist_f", frag_num_string, ".png", sep="")) 
	#dev.off()ggsave(cunt, file="~/fragmented_genome_with_snps/cunt.png")
}


