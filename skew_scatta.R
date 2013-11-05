lengths <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/skew_scatter/ex_fasta_lengths.txt", quote="\"")))
ids <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/skew_scatter/ex_ids_w_snps.txt", quote="\"")))
w <- length(lengths) #number of fragments with snps
gradients <- rep(NA, w)
abs_gradients <- rep(NA, w)
xs <- rep(NA, w) #make vectors of each x and y value
ys <- rep(NA, w)
q <- 1
for(i in lengths){
	s <- as.character(ids[q])
	p <- as.vector(as.matrix(read.table(paste("~/fragmented_genome_with_snps/skew_scatter/snps", s, ".txt", sep=""), quote="\"")))
	if (length(p) != 1){ # eliminating frags with just 1 SNP
		a <- density(p, n=(i/5)) #THE DENSITY ESTIMATE DOESNT WORK BECAUSE SOME FRAGS HAVE ONLY 1 SNP, if you only included frags with 2+ SNPs it could work
		#^ whilst n varies for each fragment, the size of each section is the same for each frag
		e <- c(which(a$x > i), which(a$x <= 0)) #vector of indices which are greater than the length or <=0. need to loop and remove each of these indices from a$x and a$y
		for (j in e){
			a$x[j] <- NA
			a$y[j] <- NA
		}
		x <- a$x[!is.na(a$x)]
		y <- a$y[!is.na(a$y)]
	}
	#xs[q] <- x
	#ys[q] <- y
	abs_gradients[q] <- abs(coef(lm(y ~ x))[2]) #gradients - NO LONGER IN SAME ORDER AS ID, ONES FROM FRAGS WITH 1 SNP MISSING
	gradients[q] <- coef(lm(y ~ x))[2]
	q <- q+1
}

png("~/fragmented_genome_with_snps/figures/figures_d2a/skew_scatter.png")
plot(ids, abs_gradients, main="The absolute gradient of SNPs for each fragment in dataset 2", xlab="Fragments", ylab="Gradient (as an absolute)") 
lines(predict(loess(abs_gradients~ids)), col="red", lwd=4)
dev.off()

png("~/fragmented_genome_with_snps/figures/figures_d2a/skew_scatter2.png")
plot(ids, gradients, main="The gradient of SNPs for each fragment in dataset 2", xlab="Fragments", ylab="Gradient")
lines(predict(loess(gradients~ids)), col="red", lwd=4)
dev.off()

#png("~/fragmented_genome_with_snps/figures/figures_d2a/example_gradient_f629.png")
#plot(xs[629], ys[629], main="Example of gradient determination for fragment 629 of dataset 2", xlab="Nucleotides", ylab="Number of SNPs")
#abline(lm(ys[629] ~ xs[629]), col="blue", lwd=4)
#dev.off()

#png("~/fragmented_genome_with_snps/figures/figures_d2a/example_gradient_f321.png")
#plot(xs[321], ys[321], main="Example of gradient determination for fragment 321 of dataset 2", xlab="Nucleotides", ylab="Number of SNPs")
#abline(lm(ys[321] ~ xs[321]), col="blue", lwd=4)
#dev.off()