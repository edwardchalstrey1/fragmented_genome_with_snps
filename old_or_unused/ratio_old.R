hm <- rnorm(35, 10000000, 5000000)
ht1a <- rnorm(1500, 5000000, 1000000)
ht1 <- ht1a[which(ht1a < 7.5e+06)] #non-recombinant region = 7.5m-12.5m
ht2a <- rnorm(1500, 15000000, 1000000)
ht2 <- ht2a[which(ht2a > 1.25e+07)] #non-recombinant region = 7.5m-12.5m
ht <- c(ht1, ht2)

het_snps <- read.table("~/fragmented_genome_with_snps/arabidopsis_datasets/ratio_dataset3/het_snps.txt", quote="\"")
hom_snps <- read.table("~/fragmented_genome_with_snps/arabidopsis_datasets/ratio_dataset3/hom_snps.txt", quote="\"")
real_ht <- as.vector(as.matrix(het_snps))
real_hm <- as.vector(as.matrix(hom_snps))

hmd <- density(hm, from=0, to=18585056)
htd <- density(ht, from=0, to=18585056)
x <- (1:512)*36298.9375
y <- hmd$y/htd$y
#plot(x, hmd$y)
#plot(x, htd$y)
#plot(x, y)

real_hmd <- density(real_hm, from=0, to=18585056)
real_htd <- density(real_ht, from=0, to=18585056)
real_y <- real_hmd$y/real_htd$y
#plot(x, real_htd$y)
#plot(x, real_hmd$y)

#plot(x, real_y, xlab="Arabidopsis chromosome 4 (nucleotides)", ylab="Ratio of Homozygous SNP Density/Heterozygous SNP Density")
#plot(x, y, xlab="Arabidopsis chromosome 4 (nucleotides)", ylab="Ratio of Homozygous SNP Density/Heterozygous SNP Density")

qqp <-qqnorm(y)
#qqp <-qqnorm(real_y)

#qqp <- qqplot(y, real_y)

#qqp <- qqplot(y, y)

corr <- cor(qqp$x,qqp$y)

library("pracma", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
findpeaks(real_y) # 