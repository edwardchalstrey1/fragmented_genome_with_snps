Fitness Methods
========================================================

Count Ratio
-------

This method counts the numbers of each type of SNP, in equally sized divisions of the genome (the R code below is repeated for homozygous and heterozygous SNPs).


```r
snp_pos <- c(1,4,5,6,9) # Example SNP positions
l <- 10 # Example genome length
div <- 2 # Example number of divisions/breaks
breaks <- c(0)
for(i in 1:div){
  breaks <- c(breaks,(l/div)*i)
  }
counts <- hist(snp_pos, breaks=breaks, plot=FALSE)$counts
counts
```

```
## [1] 3 2
```

The ratio of the homozygous/heterozygous count is then calculated for each of these divisions. A compromise is made when calculating the ratio, by adding 1 to each count. This avoids dividing by zero (or dividing zero), but results in ratios that are not accurate. The idea here is that peaks in the homozygous/heterozygous SNP density ratio, should be pronounced enough to be recognisable even when calculating the ratio in this way.

```ruby
x = 0
ratios = []
div.times do
  count_ratio = ((hm_count[x] + 1).to_f / (ht_count[x] + 1).to_f) # the measure of ratio
  ratios << count_ratio
  x+=1
end
ratios
```

"Count ratio" vectors for permutations can then be plotted against an expected count ratio vector (whilst testing the algorithm, I use the count ratio vector for the correct arrangement of contigs). Taking the absolute Pearson's correlation coefficient will give us an indication of how similar the count ratios are, and is our fitness score for this method


```r
expected_ratio <- c(1, 2, 3, 4, 5, 4, 3, 2, 1) # examples
permutation_ratio <- c(1, 3, 2, 4, 4, 5, 3, 2, 1)
score <- abs(cor(expected_ratio, permutation_ratio))
score
```

```
## [1] 0.8714
```

SNP Distance
----

The fitness ``score`` for this method, is the sum of distances between adjacent homozygous SNPs. The idea here is that the correct permutation has the lowest total distance, because the SNPs are clustered together around the causative mutation.

```ruby
homozygous_snps = [1,4,5,6,9]
score = 0
homozygous_snps.each_cons(2).map { |a,b| score+=(b-a) }
score
```

In practice, this idea is unlikely to be of use. This is because we are searching for a peak in the ratio of homozygous to heterozygous SNP density, rather than just the homozygous SNP density.

Hyp distance
----------

A "Hypothetical SNP distribution" is created that represents the ratio of homozygous to heterozygous SNPs. This is achieved by randomly assigning "SNP positions" within divisions of the genome. SNP ratios are calculated across the genome, as in the Count Ratio fitness method (see above). The number of SNPs randomly positioned within each division/break, is 10 times the calculated ratio value for that division. The idea here is that the peak in the ratio of homozygous to heterozygous SNP density, will be largest around the causative mutation, therefore the correct permutation should have the lowest total "hypothetical snp" distance.

```ruby
breaks = []
(1..ratios.length).to_a.each do |i|
	breaks << (genome_length/ratios.length.to_f)*i
end
hyp, x = [], 0
ratios.each do |ratio| 
	(ratio*10).to_i.times do
		hyp << rand(genome_length/ratios.length.to_f) + breaks[x] # random value from within the range
	end
	x+=1
end
hyp # These don't need to be unique or integers like the real SNPs, since they are just representing a distribution
```

Calling the "SNP Distance" method on the ``hyp`` array, should give a low score for well ordered permutations.

Max Density
----

The maximum density value, of a kernel density estimation for homozygous SNPs is taken as fitness score. This assumes the SNPs are clustered together around the causative mutation.


```r
homozygous_snps <- c(1,4,5,6,9) # Example
max(density(homozygous_snps)$y) # Fitness score
```

```
## [1] 0.1788
```

In practice, this idea is unlikely to be of use. This is because we are searching for a peak in the ratio of homozygous to heterozygous SNP density, rather than just the homozygous SNP density.

Max Ratio
------

The maximum value of a ratio of homozygous (kernel) density to heterozygous density, is taken as fitness score. This assumes there is a peak in the ratio of homozygous to heterozygous kernel density, when there is a peak in the actual ratio of homozygous to heterozgous SNP density.



```r
hm <- c(1,4,5,6,9) # Example of homozygous SNPs
ht <- c(1,3,5,7,9) # Example of homozygous SNPs
max(density(hm)$y/density(ht)$y) # Fitness score
```

```
## [1] 3.893
```

Max Hyp
------

This method creates the "hypothetical SNP positions" (as in Hyp Distance above), then calls the "Max Density" method on ``hyp`` instead of on the ``homozygous_snps``.
