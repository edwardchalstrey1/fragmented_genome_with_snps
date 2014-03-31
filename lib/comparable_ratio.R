#encoding: utf-8

# Input 1: Kernel density estimate of homozygous/heterozygous snp ratio
# Input 2: Dataset name string
# Input 3: Generarion that the genetic algorithm is on
# Output: Figure of the homozygous/heterozygous snp density across the genome
plot_distribution <- function(real_ratio, location, dataset, gen){
	# x <- (1:512)*36298.9375 # genome length/512 ## TODO make this method in snp_dist.rb with x argument
	x <- (1:512)*3.90625
	png(paste("~/",location,"/", dataset,"/gen_", gen, "_best_permutation_distribution.png", sep=""))
	plot(x, real_ratio, xlab="Arabidopsis chromosome 4 (nucleotides)", ylab="Ratio of Homozygous SNP Density/Heterozygous SNP Density")
	dev.off()
}