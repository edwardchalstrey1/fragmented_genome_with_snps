Genetic algorithm performance 11/3/14
========================================================

The scatter plots below show that there is little change in the [metric scores](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Comparing%20Permutations/comparing_permutations.md) for populations from random permutations (gen 0) to those selected as "fit" by the genetic algorithm. The problem here, is that the fitness method being used by the current version of [GATOC](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) does not give an accurate representation of how well ordered the contigs are. The permutation with the best fitness is "fitter" in almost every generation (closer to 0 on this scale), though oddly the average permutation seems to be worse in gen 250 than gen 125 for this [run of the algorithm](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/arabidopsis_datasets/ratio_dataset4/run5).

The fitness method being used needs to be improved.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/ratio_dataset4/run5/met_scores0-20.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/ratio_dataset4/run5/met_scores0-100.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/ratio_dataset4/run5/met_scores.png?raw=true)
