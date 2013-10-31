lengths <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/ex_fasta_lengths.txt", quote="\"")))
w <- length(lengths) #number of fragments
z <- 1:w
gradients <- rep(NA, w)
n <- 1
for(i in lengths){
	s <- as.character(z[n])
	y <- as.vector(as.matrix(read.table(paste("~/fragmented_genome_with_snps/skew_scatter/snps", s, ".txt", sep=""), quote="\"")))
	x <- 1:(length(y))
	gradients[n] <- abs(coef(lm(y ~ x))[2]) #gradient
	n <- n+1
}

png("~/fragmented_genome_with_snps/figures/skew_scatter.png")
plot(z, gradients, main="The absolute gradient of SNPs for each fragment in dataset 2", xlab="Fragments", ylab="Gradient (as an absolute)")
lo <- loess(y~x)
lines(predict(lo), col="red", lwd=4)
dev.off()

png("~/fragmented_genome_with_snps/figures/example_gradient_f629.png")
y <- as.vector(as.matrix(read.table(paste("~/fragmented_genome_with_snps/skew_scatter/snps629.txt", sep=""), quote="\"")))
x <- 1:(length(y))
plot(x, y, main="Example of gradient determination for fragment 629 of dataset 2", xlab="Nucleotides", ylab="Number of SNPs")
abline(lm(y ~ x), col="blue", lwd=4)
dev.off()