# Input: Any argument (to make function work)
# Output: Expected distribution: homozygous/heterozygous snp ratio - vector of the homozygous density estimate points divided by the heterozygous ones
comparable_ratio <- function(a){
	aaa <- a
	hm <- rnorm(35, 10000000, 5000000)
	ht1a <- rnorm(1500, 5000000, 1000000)
	ht1 <- ht1a[which(ht1a < 7.5e+06)] #non-recombinant region = 7.5m-12.5m
	ht2a <- rnorm(1500, 15000000, 1000000)
	ht2 <- ht2a[which(ht2a > 1.25e+07)] #non-recombinant region = 7.5m-12.5m
	ht <- c(ht1, ht2)
	hmd <- density(hm, from=0, to=18585056)
	htd <- density(ht, from=0, to=18585056)
	ratio <- hmd$y/htd$y
	return(ratio)
}

# Input 1: Array of heterozygous snp positions (from a fragment permutation)
# Input 2: Array of homozygous snp positions (from a fragment permutation)
# Input 3: Homozygous/heterozygous snp ratio to compare against in qq plot - vector of homozygous density estimate points divided by the heterozygous ones (from comparable_ratio)
# Output: Correlation value
get_corr <- function(het_snps, hom_snps, ratio){
	real_ht <- as.vector(as.matrix(het_snps))
	real_hm <- as.vector(as.matrix(hom_snps))
	real_hmd <- density(real_hm, from=0, to=18585056)
	real_htd <- density(real_ht, from=0, to=18585056)
	real_ratio <- real_hmd$y/real_htd$y
	qqp <- qqplot(ratio, real_ratio, plot.it=FALSE)
    #qqp <- qqnorm(real_ratio)
	return(cor(qqp$x,qqp$y))
}

# Input 1: Array of heterozygous snp positions (from a fragment permutation)
# Input 2: Array of homozygous snp positions (from a fragment permutation)
# Input 3: Homozygous/heterozygous snp ratio to compare against in qq plot - vector of homozygous density estimate points divided by the heterozygous ones (from comparable_ratio)
# Output: Distribution from data: homozygous/heterozygous snp ratio - vector of the homozygous density estimate points divided by the heterozygous ones
get_real_ratio <- function(het_snps, hom_snps, ratio){
	real_ht <- as.vector(as.matrix(het_snps))
	real_hm <- as.vector(as.matrix(hom_snps))
	real_hmd <- density(real_hm, from=0, to=18585056)
	real_htd <- density(real_ht, from=0, to=18585056)
	real_ratio <- real_hmd$y/real_htd$y
	return(real_ratio)
}

# Input 1: Distribution from data: homozygous/heterozygous snp ratio - vector of the homozygous density estimate points divided by the heterozygous ones
# Input 2: Dataset name string
# Output: Figure of the homozygous/heterozygous snp density across the genome
plot_distribution <- function(real_ratio, dataset){
	x <- (1:512)*36298.9375 # genome length/512
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/", dataset,"/best_permutation_distribution.png", sep=""))
	plot(x, real_ratio, xlab="Arabidopsis chromosome 4 (nucleotides)", ylab="Ratio of Homozygous SNP Density/Heterozygous SNP Density")
	dev.off()
}


plot_performance <- function(gen_fits, dataset){
	x <- (0:(length(gen_fits)-1))
	png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/", dataset,"/algorithm_performance.png", sep=""))
	plot(x, gen_fits, xlab="Generation", ylab="Fitness")
	dev.off()
}