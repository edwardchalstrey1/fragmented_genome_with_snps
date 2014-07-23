Data analysis for GATOC
========================================================

Here I include information on the scripts I have used to generate the analysis of my results, whilst testing my algorithm.

Distribution plot gifs
-------------

Scripts:

1. [distribution_plots.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/distribution_plots.rb)

2. [magick.sh](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/magick.sh)

Running ``ruby distribution_plots.rb dataset run gen div`` in bash shell will generate plots of the homozygous, heterozygous and ["hypothetical/ratio"](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/fitness_methods/fitness_methods.md) SNP distributions for the best permutations of each generation, in one run of the algorithm and the correct permutation. ``dataset`` is the name of the dataset directory (within ``fragmented_genome_with_snps/arabidopsis_datasets``), ``run`` is the sub directory within dataset for this run of the algorithm, ``gen`` is the number of generations you wish to make plots for, and ``div`` is the number of divisions of the genome at which to create "hypothetical SNPs".

Navigate to the run directory ``fragmented_genome_with_snps/arabidopsis_datasets/dataset/run`` and run ``bash ../../../magick.sh``. This generates gif animations for how the homozygous and "hypothetical" SNP distributions of the best permutations in each generation change as the algorithm progresses.

