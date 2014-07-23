GATOC.evolve
========================================================

Put simply, my genetic algorithm works as follows:

1. An initial population of permutations is generated, by shuffling the contig order in each. The population is an array of permutations, where each permutation is itself an array of Bio::FastaFormat objects (contigs). 

2. A fittest selection of the population is taken to be mutated (see Mutation above). A new population is created from the mutant permutations, saved fittest permutations from the previous population, and random permutations.

3. The new population undergoes selection, and another population is made. This process loops for a given number of generations, or until a quit condition is satisfied

It is used like this: ``GATOC.evolve(fasta_file, vcf_file, parameters)``. ``parameters`` is a hash of symbol keys that point to values used in the sub-methods of GATOC. Each parameter has a default value, shown in the ``opts`` hash. The parameters are described in the annotation of [GATOC.evolve](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb).

The flow chart below shows the key inputs and outputs for GATOC.evolve. For more detail on some of the sub methods see [new_population](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/new_population.md) and [select](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/select.md).

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_flowcharts/evolve.png?raw=true)