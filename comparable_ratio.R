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

qq_real_expect <- function(het_snps, hom_snps, ratio){
	real_ht <- as.vector(as.matrix(het_snps))
	real_hm <- as.vector(as.matrix(hom_snps))
	real_hmd <- density(real_hm, from=0, to=18585056)
	real_htd <- density(real_ht, from=0, to=18585056)
	real_ratio <- real_hmd$y/real_htd$y
	qqp <- qqplot(ratio, real_ratio, plot.it=FALSE)
  #qqp <- qqnorm(real_ratio)
	return(cor(qqp$x,qqp$y))
}