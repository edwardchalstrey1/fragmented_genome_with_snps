Progress 21/3/14
========================================================

Since I am in need of a new way of calculating the ratio of homozygous to heterozygous SNPs (and generating a vector that can be compared to an example distribution with a Q-Q plot), I have come up with this hypothetical SNP ratio idea:

First here I generate example list of homozygous and heterozygous SNPs (same distributions used in my model):

```r
hm <- rnorm(10000, 1e+07, 1e+06)
ht <- runif(10000, 1, 18585056)
```


Then I calculate the frequency at which SNPs occur in defined intervals:

```r
div <- 10000  # div = the number of intervals (divisions)
l <- 18585056
breaks <- c(0)
for (i in 1:div) {
    breaks <- c(breaks, (l/div) * i)
}
hmc <- hist(hm, breaks = breaks, plot = FALSE)$counts
htc <- hist(ht, breaks = breaks, plot = FALSE)$counts
```


The ratio of the frequency counts for each interval are then taken (homozygous/heterozygous). A number of "hypothetical snp positions" are added to a new array for each of the intervals. These positions are randomly chosen values within the intervals, and the number of them is dependent on the ratio frequency. The intervals are small enough, that the random choosing of points within them doesn't affect the overall distribution of hypothetical positions within the new array. The distribution of these hypothetical snps is therefore an reasonable model of the ratio distribution. These steps were done in ruby, see class [SNPdist](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/normal/lib/snp_dist.rb) methods fratio (frequency ratio) and hyp_snps (hypothetical snps).

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/normal/test/hypothetical_snps/10K_f*10/correct/d_ratio.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/normal/test/hypothetical_snps/10K_f*10/correct/hyp.png?raw=true)

What the above figures show is the density plot of the hyopthetical snp ratio is comparable to the ratio of kernel density $y values for the homozygous snp vector / the heterozygous one. This means that reforming this pattern for the hyp snps with a contig permutation, will indicate that the permutation is correct. The hyp snps plot does differ at the tails of the distribution, this however should not be a problem, as the peak is clearly identifiable. These tails are an artefact of the way I have calculated the frequency ratio, adding 1 to each of the frequencies so that a ratio can be calculated.
