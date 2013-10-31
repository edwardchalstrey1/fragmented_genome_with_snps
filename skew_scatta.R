lengths <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/ex_fasta_lengths.txt", quote="\"")))
w <- length(lengths) #number of fragments
z <- 1:w
gradients <- rep(NA, w)
n <- 1
for(i in lengths){
	s <- as.character(z[n])
	y <- as.vector(as.matrix(read.table(paste("~/fragmented_genome_with_snps/skew_scatter/snps", s, ".txt", sep=""), quote="\"")))
	x <- 1:(length(y))
	mod <- lm(y ~ x)
	gradients[n] <- abs(coef(mod)[2]) #gradient
	n <- n+1
}
png("~/fragmented_genome_with_snps/figures/skew_scatter.png")
plot(z, gradients, main="Graph to show the absolute gradient of SNPs for each fragment", xlab="Fragments", ylab="Gradient (as an absolute)")
lo <- loess(y~x)
lines(predict(lo), col="red", lwd=4)
dev.off()