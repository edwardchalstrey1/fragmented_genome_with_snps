How useful are the fitness method I am using in my genetic algorithm, and the distance metrics I am using to evaluate it's performance?
--------

Below are the results of the following experiment: I am plotting how my fitness score, and each distance metric, changes as we move away from the optimum contig order in the search space. I am taking neighbours in the search space as permutations one adjacent swap away from eachother.

To create the permutations in each population, I take the previous population, and mutate each permutation with [PMeth.adjacent_swap](https://github.com/edwardchalstrey1/pmeth). All the permutations in the starting population are identical: the optimal permutation, correctly ordered. **The size of the population I use is 10.**

In each plot, the red shaded area shows the standard deviation around the mean for each population, and the green shaded area shows the standard deviation around the mean of a population of randomly ordered permutations.

How useful are the fitness methods?
--------

### Count ratio

The plot below shows how the fitness score I use in my genetic algorithm, changes as the number of adjacent swaps from being correct a permutation is. The genome being used to generate all the plots in this document is from [small_dataset2a](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/arabidopsis_datasets/small_dataset2a), a model genome with the following parameters:


```r
# Genome length = 2Kb
# Contig lengths = 25 to 50b
# Numer of contigs = 53
# Homozyous SNP distribution = rnorm(50, 1000, 100)
# Heterozygous SNP distribution = runif(50, 1, 2000) 
```

**The number of divisions of the genome at which I calculate the ratio of homozgous to heterozygous SNPs, for the fitness score, is 100.**

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_Fitness_2000pop_10size_0.1Kdiv_swap1.png?raw=true)

What the plot above shows, is that my fitness method does consistently attribute lower scores, to permutations that are a greater number of adjacent swaps away from the correctly ordered permutation. Therefore, in my genetic algorithm, permutations being attributed with high fitness scores should be approaching the correct order. The plot shows that at 2000 adjacent swaps from the correct order, the permutations are approaching fitness scores close to that of randomly ordered permutations, but are not there yet, and the rate of decrease appears to be slowing.

### SNP distance

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/thirdfitness/arabidopsis_datasets/small_dataset2b/adjacent_swaps_Fitness_2000pop_10size_swap1_snp_distance.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/thirdfitness/arabidopsis_datasets/small_dataset2b/adjacent_swaps_Fitness_10000pop_10size_swap1_snp_distance.png?raw=true)

### Max density (homozygous SNPs)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/thirdfitness/arabidopsis_datasets/small_dataset2b/adjacent_swaps_Fitness_2000pop_10size_swap1_max_density.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/max_density/arabidopsis_datasets/small_dataset2b/adjacent_swaps_Fitness_10000pop_10size_swap1_density.png?raw=true)

### Max ratio (homozygous/heterozygous SNPs)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/thirdfitness/arabidopsis_datasets/small_dataset2b/adjacent_swaps_Fitness_2000pop_10size_swap1_max_ratio.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/thirdfitness/arabidopsis_datasets/small_dataset2b/adjacent_swaps_Fitness_10000pop_10size_swap1_max_ratio.png?raw=true)


### Max hyp (hypothetical distribution based on ratio of homozgous to heterozygous SNPs)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2c/adjacent_swaps_Fitness_2000pop_10size_swap1_max_hyp.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2c/adjacent_swaps_Fitness_10000pop_10size_swap1_max_hyp.png?raw=true)

### Hyp distance

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2c/adjacent_swaps_HypDistance_2000pop_10size_swap1.png)

How useful are distance metrics
--------

### Deviation distance

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_DeviationDistance_2000pop_10size_0.1Kdiv_swap1.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_DeviationDistance_10000pop_10size_0.1Kdiv_swap1.png?raw=true)

### Squared deviation distance

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_SquareDeviationDistance_2000pop_10size_0.1Kdiv_swap1.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_SquareDeviationDistance_10000pop_10size_0.1Kdiv_swap1.png?raw=true)

### Generalized hamming distance

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_HammingDistance_2000pop_10size_0.1Kdiv_swap1.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_HammingDistance_10000pop_10size_0.1Kdiv_swap1.png?raw=true)

### R distance

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_RDistance_2000pop_10size_0.1Kdiv_swap1.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_RDistance_10000pop_10size_0.1Kdiv_swap1.png?raw=true)

### Longest common sub-sequence

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_LongestCommonSubsequence_2000pop_10size_0.1Kdiv_swap1.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_LongestCommonSubsequence_10000pop_10size_0.1Kdiv_swap1.png?raw=true)

### Kendall's Tau

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_KendallsTau_2000pop_10size_0.1Kdiv_swap1.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset2a/adjacent_swaps_KendallsTau_10000pop_10size_0.1Kdiv_swap1.png?raw=true)
