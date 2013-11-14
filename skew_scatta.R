lengths <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/skew_scatter/ex_fasta_lengths.txt", quote="\""))) #frags_w_snps
ids <- as.vector(as.matrix(read.table("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/skew_scatter/ex_ids_w_snps.txt", quote="\""))) #frags_w_snps
w <- length(lengths) #number of fragments with snps
gradients <- c()
abs_gradients <- c()
xs <- c() #make vectors of each x and y vector
ys <- c()
q <- 1
n <- 1
nu_ids <- c()
for(i in lengths){ #for each fragment
	s <- as.character(ids[q])
	p <- as.vector(as.matrix(read.table(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/skew_scatter/snps", s, ".txt", sep=""), quote="\"")))
	if (length(p) > 6){ # eliminating frags with less than 6 SNPs
		nu_ids <- c(nu_ids, ids[q]) # add the id to a new vector, if it has 6 or more SNPs
		a <- density(p, n=i/5) #dividing the length of the fragment by a constant (5) to get the number of equally spaced points at which the density is to be estimated
		#^ whilst n varies for each fragment, the size of each section is the same for each frag
		e <- c(which(a$x > i), which(a$x <= 0)) #vector of indices which are greater than the length or <=0. need to loop and remove each of these indices from a$x and a$y
		for (j in e){
			a$x[j] <- NA # making those redundant indices NA
			a$y[j] <- NA
		}
		x <- a$x[!is.na(a$x)] #removing indices that are NA
		y <- a$y[!is.na(a$y)]

		xs[[n]] <- x
		ys[[n]] <- y
		n <- n+1

		coeff <- coef(lm(y ~ x))
		abs_gradients <- c(abs_gradients, abs(coeff[2])) #gradients - should now be in same order as nu_ids
		gradients <- c(gradients, coeff[2])
	}
	q <- q+1
}


png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/skew_scatter.png")
plot(nu_ids, abs_gradients, main="The absolute gradient of SNPs for fragments with 6 or more SNPs from dataset 6", xlab="Fragments", ylab="Gradient (as an absolute)") 
lines(lowess(nu_ids, abs_gradients), col="red", lwd=4)
dev.off()

png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/skew_scatter2.png")
plot(nu_ids, gradients, main="The gradient of SNPs for fragments with 6 or more SNPs from dataset 6", xlab="Fragments", ylab="Gradient")
lines(lowess(nu_ids, gradients), col="red", lwd=4)
dev.off()

png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/example_gradient_f687.png")
z <- which(687 == nu_ids)
plot(xs[[z]], ys[[z]], main="Example of gradient determination for fragment 687 of dataset 6", xlab="Fragment", ylab="SNP density (Kernel density estimation)") 
dev.off()

png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/example_gradient_f257.png")
zz <- which(257 == nu_ids)
plot(xs[[zz]], ys[[zz]], main="Example of gradient determination for fragment 257 of dataset 6", xlab="Fragment", ylab="SNP density (Kernel density estimation)")
dev.off()

png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/example_gradient_f1042.png")
zzz <- which(1042 == nu_ids)
plot(xs[[zzz]], ys[[zzz]], main="Example of gradient determination for fragment 1042 of dataset 6", xlab="Fragment", ylab="SNP density (Kernel density estimation)")
dev.off()


snp_pos <- as.vector(as.matrix(read.table(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/skew_scatter/snp_pos.txt", sep=""), quote="\"")))
snp_dens <- density(snp_pos, n=(sum(lengths)/5))
e <- c(which(snp_dens$x > sum(lengths)), which(snp_dens$x <= 0)) #vector of indices which are greater than the genome length or <=0.
for (j in e){
	snp_dens$x[j] <- NA # making those redundant indices NA
	snp_dens$y[j] <- NA
}
x <- snp_dens$x[!is.na(snp_dens$x)] #removing indices that are NA
y <- snp_dens$y[!is.na(snp_dens$y)]

png("~/fragmented_genome_with_snps/arabidopsis_datasets/dataset6/figures/genome_kdens_snps.png")
plot(x, y, main="The Kernel density estimation for SNPs across the model genome (dataset 6)", xlab="Arabidopsis chromosome 4 (bases)", ylab="SNP density (Kernel density estimation)")
dev.off()