How useful is the fitness method for finding permutations close to the correct order
========================================================

The plot below shows how the fitness score I use in my genetic algorithm, changes as we move away from the optimum contig order in the search space. The genome being used to generate all the plots in this document is from [small_dataset2a](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/arabidopsis_datasets/small_dataset2a), a model genome with the following parameters:


```r
# Genome length = 2Kb
# Contig lengths = 25 to 50b
# Numer of contigs = 53
# Homozyous SNP distribution = rnorm(50, 1000, 100)
# Heterozygous SNP distribution = runif(50, 1, 2000) 
```

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_Fitness_100.png?raw=true)


How useful are distance metrics
======

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_DeviationDistance_100pop_10size_0.1Kdiv.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_SquareDeviationDistance_100pop_10size_0.1Kdiv.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_HammingDistance_100pop_10size_0.1Kdiv.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_RDistance_100pop_10size_0.1Kdiv.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_LongestCommonSubsequence_100pop_10size_0.1Kdiv.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_KendallsTau_100pop_10size_0.1Kdiv_not_blank.png?raw=true)
