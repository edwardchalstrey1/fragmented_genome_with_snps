Change in homozygous to heterozygous SNP ratio over generations of genetic algorithm
========================================================

To create the gifs:
-------

1. Create the plots with ```ruby distribution_plots.rb```
 - dataset = ARGV[0]
 - run = ARGV[1]
 - gen = ARGV[2]
 - div = ARGV[3] # number of divisions for the genome that SNP ratios calculated in this run

2. From inside the run directory: ```bash ../../../magick.sh```

How does the distribution change for the small dataset
---------

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset4/run1/images_hyp.gif?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset4/run1/Gencorrect_lists/best_permutation_distribution_hyp_0.1Kdiv.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset4/run1/images_ratios_0.1Kdiv.gif?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/small_dataset4/run1/Gencorrect_lists/best_permutation_ratios_0.1Kdiv.png?raw=true)

How does the distribution change for the main dataset (10K_dataset4) parallel run 60
------------

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run60/images_hyp.gif?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run60/Gencorrect_lists/best_permutation_distribution_hyp_10.0Kdiv.png?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run60/images_ratios_10Kdiv.gif?raw=true)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run60/Gencorrect_lists/best_permutation_ratios_10Kdiv.png?raw=true)

How does the distribution change for the main dataset (10K_dataset4) parallel run 12
--------

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run12/images_hyp.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run12/Gencorrect_lists/best_permutation_distribution_hyp_1.0Kdiv.png?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run12/images_ratios.gif?raw=true)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/10K_dataset4/p_run12/Gencorrect_lists/best_permutation_distribution_hyp_1.0Kdiv.png?raw=true)
