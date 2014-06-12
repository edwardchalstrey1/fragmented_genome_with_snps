Detecting Causative Mutations from High-throughput Sequencing on Unordered Genomes
===========================

[![Build Status](https://drone.io/github.com/edwardchalstrey1/fragmented_genome_with_snps/status.png)](https://drone.io/github.com/edwardchalstrey1/fragmented_genome_with_snps/latest)

I am creating a model genome, based of [Arabidopsis chromosome 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/TAIR10_chr4.fasta), from an individual created by the back-crossing of a mutant with a mapping line, where the mutation is an experimentally induced, phenotype altering SNP. 

I am designing an algorithm that can locate a causative SNP mutation, based on the distribution of SNPs. To emulate real data, I first split the genome into fragments, which model contigs assembled from high-throughput sequencing reads. My algorithm will need to rearrange these contigs into their proper order, then locate the causative mutation.

### The Experiment Being Modelled
 
A beneficial Arabidopsis phenotype is created by experimentally by inducing SNPs with EMS (e.g. a disease resistance mutation). The phenotype altering, reccesive homozygous SNP needs to be identified. The mutant individual is back-crossed with a mapping line. This mapping line is a cross between two Arabidopsis ecotypes (one of which the same ecotype as the mutant), and therefore has heterozygous SNPs across its genome.

Back-crossing like this multiple times, but always selecting for progeny with the mutant phenotype, will result in an individual that has a non-recombinant region with high homozygous to heterozygous SNP ratio near the mutant position, due to linkage. This is because with each back-cross, the EMS induced homozygous SNPs closest to mutation being selected for have a higher chance of being conserved through linkage, than those further from the location of the mutation. 

In other words, the remaining linked homozygous SNPs will be distributed around the phenotype causing SNP. Plotting the homozygous/heterozygous SNP ratio, should reveal a peak in the non-recombinant region, because of the high density of homozygous SNPs in that location (there is a consistent distribution of heterozygous SNPs from the mapping line across the genome). The causative SNP should be located at the peak of this ratio distribution.

### Creating a model genome

Running: **ruby [create_model_genome.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/create_model_genome.rb) ARGV[0] = dataset_name** will generate a new model genome of that name based on Arabidopsis chromosome 4 and the experiment detailed above, in fragmented_genome_with_snps/arabidopsis_datasets. This includes a FASTA file with the sequences of each fragment, and a VCF file with the SNPs on each fragment. In the INFO field of the VCF, each SNP has been given an allele frequency (AF). Heterozygous SNPs will generally have AF = ~0.5, and homozygous AF = ~1.0, but this will vary with pooled data. In the model, each SNP has been given an allele frequency of exactly 0.5 or 1.0. The methods used to create these files are in the [ModelGenome class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/arabidopsis_c4_w_snps.rb).

Why not use real data straight away?
 - With model data, I know the location of the causative mutation and the correct order of the fragments. This means I can tell whether or not my algorithm works.
 - I can generate datasets that model different experiments/sample data quickly: e.g. different fragment sizes, different SNP distributions, different experiment

Re-ordering the Genome
----------

I aim to use a **genetic algorithm** to rearrange contigs, and find a permutation of the contig order that has the correct (expected) homozygous/heterozygous SNP ratio distribution. Permutations with distributions close to the expected, are likely be ordered in an arrangement close to the correct. The algorithm I am creating will eventually be capable of determining the location of the causative mutation, in a genome from the experiment detailed above. The algorithm will take a FASTA file of all the unordered genome fragments, and a VCF file containing the SNP positions for each fragment.

### Genetic Algorithm

A genetic algorithm is a kind of iterative improvement algorithm. These are where an initial starting solution to a problem is incrementally improved, in order to reach an eventual perfect solution. In the case of a permutation problem, like my contig rearrangement, it is the *ordering* of elements (contigs) that needs to be changed until correct.

Genetic algorithms are based on the principles of natural selection. A "population" of solutions is created, and those with the highest "fitness" (correctness) are selected. New "offspring" solutions can be created through recombination of the components of "parent" solutions, and "mutations" can be introduced by creating novel components or parameters. If subsequent populations are created by the recombination of selected members of the current population, and/or the introduction of mutant solutions, the best solutions in each generation will be better than the best solutions from the previous generation. In this way, the perfect solution (in this case permutaion of fragments) will eventually be reached.

With my model genome, I know the correct permutation of contigs/fragments (I have named each frag1, frag2...), so I am able to see if the algorithm has arranged them correctly (see 'Assesing the algorithms performance' below). I give permutations of the fragment order a fitness score based on how well the expected homozygous/heterozygous SNP ratio distribution is reformed (see permutation fitness below).

### Mutation

In my genetic algorithm, I use mutation methods to rearrange the contigs into novel permutations. The mutation methods used are from my [PMeth ruby gem](https://github.com/edwardchalstrey1/pmeth). The repo provides a description of what they do.
In each generation of my genetic algorithm, I create a number of mutants from the previous population, but also save some the permutations with best fitness scores. I also include a number of random permutations in each generation.

### Permutation Fitness

Genetic algorithms require a fitness method, in order to tell how close each solution comes to solving the given problem. With a permutation problem, the fitness method should score the permutation based on how correct the ordering of elements is. In my algorithm, I want to identify a mutation based on the ratio of homozygous to heterozygous SNP distributions in a genome. I will therefore use the ratio as a way of telling how close a given fragment permutation is to the correct order. My algorithm does this by working out where the SNP positions are in the genome, and assuming the contigs are ordered in the way of a given permutation (the position of SNPs on each fragment are known). The idea here, is that the closer the permutation is to being correct, the closer the ratio distribution will be to the expected distribution.

The ``score`` method in [FitnessScore](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/snp_dist.rb) is used within the fitness method of my genetic algorithm. It compares ratio of homozygous to heterozygous SNP counts of a contig permutation, with that of the correct permutation. It does this by plotting the vectors of counts (for the permutation and correct arrangement) against eachother. The Pearson correlation coefficient plot is obtained; a value between 0 and 1. The closer the value is to 1, the more closer the permutation is to the correct fragment arrangement.

### Evolving the Population

The ``evolve`` method in the [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) uses all the other methods and works as follows:

1. An initial population of permutations is generated, by shuffling the contig order in each (with ruby's .shuffle method for arrays). The population is an array of permutations, where each permutation is itself an array of Bio::FastaFormat objects (contigs). 

2. A fittest selection of the population is taken to be mutated (see Mutation above). A new population is created from the mutant permutations, saved fittest permutations from the previous population, and random permutations.

3. The new population undergoes selection, and another population is made. This process loops for a given number of generations, or until the [quit condition](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Quit_condition/quit_condition.md) is satisfied.

### Finding the Peak: The Causative Mutation

I am yet to implement this into my algorithm as of 02/06/14

One possibility: Using the R library pracma, the findpeaks function can be used to identify the peaks in a distribution vector.

Running the algorithm with current project files: 02/06/14
--------

1. To run it: **ruby [algorithm_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_test.rb) ARGV[0] = dataset name ARGV[1] = run**, where the dataset name is the same as the one used to create the model genome, and the run is the name of the directory to save performance figures to. **Subsequent required command line parameters are described in [algorithm_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_test.rb).**

1. The [ReformRatio class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) contains methods involved in interpreting data from FASTA and VCF files, so they can be used in the genetic algorithm.

2. The [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) (Genetic Agorithm To Order Contigs) contains all the genetic algorithm methods. See below for genetic algorithm description.
 - Most of the other lib files contain methods used by methods within GATOC, or are used to plot figures

Assesing the Genetic Algorithm's performance
---------------------------------------

My experiment relies on the premise that the fitness of a permutation, as I have calculated it, is related to how close that permutation is to the correct arrangement of contigs. To test this, I use a number of distance metrics, which compare the order of permutations and provide a similarity score. I have made a ruby gem called [PDist](https://github.com/edwardchalstrey1/pdist) that contains the methods for calculating these scores, and the repo describes how they work.

Project Background Information
------------

Identifying the genes and alleles associated with beneficial and deleterious traits, is of the utmost importance to agronomic and biomedical science. When beneficial mutations are identified in agronomically important plant species, breeders can select for progeny with the mutant phenotype. This molecular breeding strategy speeds up the selection process, because rejected plants need not be grown. If breeders wish to enhance a complex trait, such as pathogenic disease resistance in a crop species, natural variation can be a limiting factor. One approach is to induce mutations experimentally, using mutagens such as ethylmethane sulphonate (EMS).

Identifying causal mutants is important in understanding the mechanism by which they confer a beneficial phenotype. Software packages like SHOREmap (Schneeberger et al, 2009) and NGM (Austin et al, 2011) have been developed, which can be used to identify causative mutations in F2 EMS mutants. NGM for example, works using a chastity statistic to quantify the relative contribution of the parental mutant, and mapping lines, to each single nucleotide polymorphism (SNP), in the pooled F2 population.

However, in species where a reference genome is not available, there exists a need to develop tools which can identify mutations directly from HTS datasets.

Project dependencies
------------

1. Ruby >= 2.0.0
2. Ruby gems:
 - pmeth >= 1.0.0
 - pdist >= 0.0.5
 - bio >= 1.4.3.0001
 - bio-samtools >= 2.0.5
 - rinruby >= 2.0.3
3. R >= 3.0.1
4. R packages:
 - pracma >= 1.6.4
 - ggplot2 >= 0.9.3.1
