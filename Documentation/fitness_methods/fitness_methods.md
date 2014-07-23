Fitness Methods
========================================================

Count ratio
-------


```r
correct_count_ratio <- c(1, 2, 3, 4, 5, 4, 3, 2, 1) # examples
permutation_count_ratio <- c(1, 3, 2, 4, 4, 5, 3, 2, 1)
abs(cor(correct_count_ratio, permutation_count_ratio))
```

```
## [1] 0.8714
```

SNP distance
----

```ruby
homozygous_snps = [1,4,5,6,9]
score = 0
homozygous_snps.each_cons(2).map { |a,b| score+=(b-a) }
puts score
```

The sum of the distances between the homozygous SNPs. The idea here is that the correct permutation has the lowest total distance, because the SNPs are clustered together around the causative mutation.

Max density
----


```r
homozygous_snps <- c(1,4,5,6,9) # Example
max(density(homozygous_snps)$y) # Fitness score
```

```
## [1] 0.1788
```

Max Ratio
------


```r
hm <- c(1,4,5,6,9) # Example of homozygous SNPs
ht <- c(1,3,5,7,9) # Example of homozygous SNPs
max(density(hm)$y/density(ht)$y) # Fitness score
```

```
## [1] 3.893
```
