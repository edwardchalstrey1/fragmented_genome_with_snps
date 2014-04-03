Progress 3/4/14
========================================================

I have decided to create a [smaller test dataset](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/test/arabidopsis_datasets/small_dataset2) to run my genetic algorithm on, to test whether it will work.
This is a model genome of length **2Kb** (much smaller than my main arabidopsis chr4 length ≈ 18Mb), that has been fragmented into contigs of **25 - 50b**. The SNP distributions used are as follows:


```r
hm <- rnorm(50, 1000, 100)
ht <- runif(50, 1, 2000)
```


Here I have continued to use the [hypothetical SNPs idea](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/test/Progress/Hypothetical_SNP_ratio/Hypothetical_SNP_ratio.md) as a comparable ratio for the Q-Q plot in the fitness method. Firstly I use the same SNP positions used to create the model genome, to create a list of hypothetical SNPs that match the ratio distribution (homozygous/heterozygous SNPs) across the correctly ordered genome:

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/test/arabidopsis_datasets/small_dataset2/gen_correct_ordered_contigs_best_permutation_distribution.png?raw=true)

This plot shows the ratio approximated by the hypothetical SNPs, for lists of homozygous and heterozygous SNPs derived from correctly ordered contigs. The algorithm appears to reform this ratio incredibly well, pretend that:

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/test/arabidopsis_datasets/small_dataset2/run8/gen_12_best_permutation_distribution.png?raw=true)

This goes to show that the genetic algorithm method can find permutations that have a matching SNP ratio to an example ratio. The best permutations outputted by the algorithm however, are not correctly ordered. The figure below shows that the permutation used to create the figure above had a fitness comparable to the correct permutation, but a poor square deviation distance (the same for other permutation metrics).

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/test/arabidopsis_datasets/small_dataset2/run8/square_deviation_distance_gen_0-10.png?raw=true)

The good news is that this dataset is very small, which means there are likely many permutations that place the SNPs similarly in the 2Kb genome. This is less likely to be the case for a more complex genome with a larger number of SNPs, i.e. there is a smaller chance of dramatically differently ordered contigs reforming the ratio distribution, like has happened here. This test dataset does show that in principle, the algorithm can be used to match distributions of the same shape.

**So far I have not run the most recent version of the genetic algorithm on the Arabidopsis chr4 model I have developed to completion, and will need to do so in the cluster.**

