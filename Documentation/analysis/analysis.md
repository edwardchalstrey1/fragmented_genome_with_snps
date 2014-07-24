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

Do distance metrics and fitness methods work?
-----------

See [Progress/do_metrics_work](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Do_metrics_work/do_metrics_work.md) for results and details.

Script: [do_metrics_work.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/do_metrics_work.rb)

Running ``ruby do_metrics_work.rb dataset size pop_num data_plot metric swap_num div`` will generate a plot of the change in the metric, as the distance from the correct permutations increases (by adjacent swaps of contigs). ``dataset`` is the same as above. ``size`` is the number of permutations to include in each population. ``pop_num`` is the number of populations to create. ``data_plot`` should be either 'csv', 'plot', or both. ``metric`` should be the name of the distance metric to test (from [PDist](https://github.com/edwardchalstrey1/pdist)), or the fitness method to test (from [FitnessScore](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/fitness_score.rb)), but should be camel-cased without underscores/spaces, e.g. 'kendalls_tau' metric is used when ``metric`` = 'KendallsTau'. ``swap_num`` is the number of swap mutations to be carried out on permutations between consecutive populations (I have always set this to 1). ``div`` is the number of divisions of the genome at which to calculate the ratio when testing the 'count_ratio', 'max_hyp' and 'hyp_distance' fitness methods. ``div`` can be left blank for all other methods.

