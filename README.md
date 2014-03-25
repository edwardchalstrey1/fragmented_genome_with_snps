Detecting Causative Mutations from High-throughput Sequencing on Unordered Genomes
===========================

I am creating a model genome, based of [Arabidopsis chromosome 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/TAIR10_chr4.fasta), from an individual created by the back-crossing of a mutant with a mapping line, where the mutation is an experimentally induced, phenotype altering SNP. 

I am designing an algorithm that can locate a causative SNP mutation, based on the distribution of SNPs. To emulate real data, I first split the genome into fragments, which model contigs assembled from high-throughput sequencing reads. My algorithm will need to rearrange these fragments into their proper order, then locate the causative mutation.

Creating a model genome
--------------------

Running: **ruby [create_model_genome.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/create_model_genome.rb) ARGV[0] = dataset_name** will generate a new model genome of that name based on Arabidopsis chromosome 4 and the experiment detailed above, in fragmented_genome_with_snps/arabidopsis_datasets. This includes a FASTA file with the sequences of each fragment, and a VCF file with the SNPs on each fragment. In the INFO field of the VCF, each SNP has been given an allele frequency (AF). Heterozygous SNPs will generally have AF = ~0.5, and homozygous AF = ~1.0, but this will vary with pooled data. In the model, each SNP has been given an allele frequency of exactly 0.5 or 1.0. The methods used to create these files are in the [ModelGenome class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/arabidopsis_c4_w_snps.rb).

Why not use real data straight away?
 - With model data, I know the location of the causative mutation and the correct order of the fragments. This means I can tell whether or not my algorithm works.
 - I can generate datasets that model different experiments/sample data quickly: e.g. different fragment sizes, different SNP distributions, different experiment
 
### The Experiment Being Modelled
 
A beneficial Arabidopsis phenotype is created by experimentally by inducing SNPs with EMS (e.g. a disease resistance mutation). The phenotype altering, reccesive homozygous SNP needs to be identified (cannot be hetetozygous or the selection experiment wouldn't work). The mutant individual is back-crossed with a mapping line. This mapping line is a cross between two Arabidopsis ecotypes (one of which the same ecotype as the mutant), and therefore has heterozygous SNPs across its genome.

Back-crossing like this multiple times, but always selecting for progeny with the mutant phenotype, will result in an individual that has a non-recombinant region with high homozygous to heterozygous SNP ratio near the mutant position, due to linkage. This is because with each back-cross, the EMS induced homozygous SNPs closest to mutation being selected for have a higher chance of being conserved through linkage, than those further from the location of the mutation. 

In other words, the remaining linked homozygous SNPs will be (normally) distributed around the phenotype causing SNP. Plotting the homozygous/heterozygous SNP ratio should reveal the non-recombinant region has a high ratio, because of the high density of homozygous SNPs here, when compared to the consistent distribution of heterozygous SNPs (from the mapping line); the causative SNP should be located at the peak of this ratio distribution.

Re-ordering the Genome
----------

I am currently attempting to use a **genetic algorithm** to rearrange the fragments, by reforming the homozygous/heterozygous SNP ratio distribution. The algorithm I am creating will eventually be capable of determining the location of the causative mutation, in a genome from the experiment detailed above. The algorithm will take a FASTA file of all the unordered genome fragments, and a VCF file containing the SNP positions for each fragment. Currently, the algorithms below can only be used for my model data, based on Arabidopsis c4 (see above).

1. To run it: **ruby [algorithm_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_test.rb) ARGV[0] = dataset name ARGV[1] = run **, where the dataset name is the same as the one used to create the model genome, and the run is the name of the directory to save performance figures to.

1. The [ReformRatio class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) contains methods involved in interpreting data from FASTA and VCF files, so they can be used in the genetic algorithm.

2. The [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) (Genetic Agorithm To Order Contigs) contains all the genetic algorithm methods. See below for genetic algorithm description.

3. [comparable_ratio.R](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/comparable_ratio.R) is used within the fitness method of my genetic algorithm. It compares the homozygous/heterozygous SNP distribution of the rearranged fragments, to the same distributions used when creating the model genome, using a qq plot. The Pearson correlation coefficient (r) for the Q-Q plot is obtained, between 0 and 1. The closer the value is to 1, the more likely it is that the correct fragment arrangement has been found.

Genetic Algorithm
-------

Using a genetic algorithm, I will find the correct permutation of fragments to reform the genome. With my model genome, I know the correct permutation of fragments (I have named each frag1, frag2...), so I am able to see if the algorithm has arranged them correctly. I will score permutations of the fragment order by how well the expected homozygous/heterozygous SNP ratio distribution is reformed (see permutation fitness below).

A genetic algorithm is a kind of iterative improvement algorithm. These are where an initial starting solution to a problem is incrementally improved, in order to reach an eventual perfect solution. In the case of a permutation problem, like my fragment rearrangement, it is the ordering of elements (fragments) that needs to be changed until correct.

Genetic algorithms are based on the principles of natural selection. A "population" of solutions is created, and those with the highest "fitness" (correctness) are selected. New "offspring" solutions can be created through recombination of the components of "parent" solutions, and "mutations" can be introduced by creating novel components or parameters. If subsequent populations are created by the recombination of selected members of the current population, and the introduction of mutant solutions, the idea is that the best solutions in each generation will be better than the best solutions in the previous generation. In this way, the perfect solution (in this case permutaion of fragments) will eventually be reached.

### Recombination

The recombination method of the genetic algorithm in the [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) works as follows:
Two parent permutations are used as input. Each of the parents are split into chunks of size X, where X is a random whole number that the number of fragments is divisible by. One of the chunks from parent 2 is chosen at random, to be recombined with parent 1. The "child" permutation, starts out as a copy of parent 1, then the chosen chunk from parent 2 replaces the equivalent chunk from parent 1. To avoid duplicating fragments (and losing others), each fragment from the chunk being displaced by the "chosen chunk", is placed into the position it's corresponding fragment (that holds the same position in the chosen chunk) occupies in parent 1. This results in a child permutation that has some parts ordered in the same way as parent 1, a chunk that is in the same order as in parent 2, and some fragments different positions than in either parent.

### Mutation

The mutation method of the genetic algorithm in the [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) works as follows:
One permutation is taken as input (the permutation to be mutated). It is split into chunks of size X, where X is a random whole number that the number of fragments is divisible by (just like in the recombination method). The fragments within the chunk are then shuffled, to give the new "mutant" permutation.

### Permutation Fitness: Q-Q plot

Genetic algorithms require a fitness method, in order to tell how close each solution comes to solving the given problem. With a permutation problem, the fitness method should score the permutation based on how correct the ordering of elements is. In my algorithm, I want to identify a mutation based on the ratio of homozygous to heterozygous SNP distributions in a genome. I will therefore use the ratio as a way of telling how close a given fragment permutation is to the correct order. My algorithm does this by working out where the SNP positions are in the genome, assuming the fragments are ordered in the way of a given permutation (the position of SNPs on each fragment are known). The idea here, is that the more correct the arrangement of fragments, the closer the ratio distribution will be to the expected distribution.

To calculate how similar the reformed SNP ratio (from a given fragment permutation) is to the expected ratio, a Q-Q plot is used. A Q-Q plot is a graphical method for comparing two probability distributions by plotting their quantiles against each other. By plotting the reformed distribution against the expected distribution, a correleation value is obtained (see Re-ordering the genome: [comparable_ratio.R](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/comparable_ratio.R) above). This is the value returned by the fitness method.

### Evolving the Population

The evolve method in the [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) works as follows:
First an initial population of fragment order permutations is generated, by shuffling the fragments. The number of permutations in a population (the population size) is user defined. Then fittest half (or user defined amount) is selected for recombination. Each of the selected permutations is recombined with another randomly chosen permutation (out of the selection), and a user defined number of mutants are created (currently based on the fittest permutation from the previous generation). A user defined number of the fittest permutations are also carried over into the next generation. What this means is that the fittest permutation in each generation will either be more fit, or the same as the fittest permutation in the previous generation.

### Finding the Peak: The Causative Mutation

USE PEAK FINDING IN R

Using the R library pracma, the findpeaks function can be used to identify the peaks in a distribution vector.

Assesing the Genetic Algorithm's performance
---------------------------------------

To determine how close each permutation outputted by the genetic algorithm is to being correct, I use a number of [distance metrics](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/normal/Comparing%20Permutations/comparing_permutations.md). To generate figures of how the algorithm performs over generations with regard to fitness and these metrics run: **ruby [ga_performance.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/normal/ga_performance.rb) ARGV[0] = dataset name ARGV[1] = run **, where the dataset name is the same as the one used to create the model genome, and the run is the run of the algorithm you wish to generate the figures for.

Project Background Information
------------

Identifying the genes and alleles associated with beneficial and deleterious traits, is of the utmost importance to agronomic and biomedical science. When beneficial mutations are identified in agronomically important plant species, breeders can select for progeny with the mutant phenotype. This molecular breeding strategy speeds up the selection process, because rejected plants need not be grown. If breeders wish to enhance a complex trait, such as pathogenic disease resistance in a crop species, natural variation can be a limiting factor. One approach is to induce mutations experimentally, using mutagens such as ethylmethane sulphonate (EMS).

Identifying causal mutants is important in understanding the mechanism by which they confer a beneficial phenotype. Software packages like SHOREmap (Schneeberger et al, 2009) and NGM (Austin et al, 2011) have been developed, which can be used to identify causative mutations in F2 EMS mutants. NGM uses a chastity statistic to quantify the relative contribution of the parental mutant, and mapping lines, to each single nucleotide polymorphism (SNP), in the pooled F2 population.

However, in species where a reference genome is not available, there exists a need to develop tools which can identify mutations directly from HTS datasets.


