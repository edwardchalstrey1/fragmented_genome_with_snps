Alternate fitness method tests with small dataset
========================================================

SNP distance
----

```ruby
score = 0
homozygous_snps.each_cons(2).map { |a,b| score+=(b-a) }
```

The sum of the distances between the homozygous SNPs. The idea here is that the correct permutation has the lowest total distance, because the SNPs are clustered together around the causative mutation.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/p_run1/images_hm.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/p_run1/Gencorrect_lists/best_permutation_distribution_hm_0.1Kdiv.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/p_run1/images_hyp.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/p_run1/Gencorrect_lists/best_permutation_distribution_hyp_0.1Kdiv.png?raw=true)

Max density
----


```r
homozygous_snps <- c(1,4,5,6,9) # Example
max(density(homozygous_snps)$y) # Fitness score
```

```
## [1] 0.1788
```

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxd2/images_hm.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxd2/Gencorrect_lists/best_permutation_distribution_hm_0.1Kdiv.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxd2/images_hyp.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxd2/Gencorrect_lists/best_permutation_distribution_hyp_0.1Kdiv.png?raw=true)

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

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxr/images_hm.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxr/Gencorrect_lists/best_permutation_distribution_hm_0.1Kdiv.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxr/images_hyp.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/testmaxr/Gencorrect_lists/best_permutation_distribution_hyp_0.1Kdiv.png?raw=true)
