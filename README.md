GATOC (Genetic Algorithm To Order Contigs)
===========================

Background Information
------------

Identifying the genes and alleles associated with beneficial and deleterious traits, is of the utmost importance to agronomic and biomedical science. When beneficial mutations are identified in agronomically important plant species, breeders can select for progeny with the mutant phenotype. This molecular breeding strategy speeds up the selection process, because rejected plants need not be grown. If breeders wish to enhance a complex trait, such as pathogenic disease resistance in a crop species, natural variation can be a limiting factor. One approach is to induce mutations experimentally, using mutagens such as ethylmethane sulphonate (EMS).

Identifying causal mutants is important in understanding the mechanism by which they confer a beneficial phenotype. Software packages like SHOREmap (Schneeberger et al, 2009) and NGM (Austin et al, 2011) have been developed, which can be used to identify causative mutations in F2 EMS mutants. NGM for example, works using a chastity statistic to quantify the relative contribution of the parental mutant, and mapping lines, to each single nucleotide polymorphism (SNP), in the pooled F2 population.

However, in species where a reference genome is not available, there exists a need to develop tools which can identify mutations directly from HTS datasets.

Project
-----

I am creating a model genome, based of [Arabidopsis chromosome 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/TAIR10_chr4.fasta), from an individual created by the back-crossing of a mutant with a mapping line, where the mutation is an experimentally induced, phenotype altering SNP. 

I am designing an algorithm that can locate a causative SNP mutation, based on the distribution of SNPs. To emulate real data, I first split the genome into fragments, which model contigs assembled from high-throughput sequencing reads. My algorithm will need to rearrange these contigs into their proper order, to locate the causative mutation.

### The Experiment Being Modelled
 
A beneficial Arabidopsis phenotype is created by experimentally by inducing SNPs with EMS (e.g. a disease resistance mutation). The phenotype altering, reccesive homozygous SNP needs to be identified. The mutant individual is back-crossed with a mapping line. This mapping line is a cross between two Arabidopsis ecotypes (one of which the same ecotype as the mutant), and therefore has heterozygous SNPs across its genome.

Back-crossing like this multiple times, but always selecting for progeny with the mutant phenotype, will result in an individual that has a non-recombinant region with high homozygous to heterozygous SNP ratio near the mutant position, due to linkage. This is because with each back-cross, the EMS induced homozygous SNPs closest to mutation being selected for have a higher chance of being conserved through linkage, than those further from the location of the mutation. 

In other words, the remaining linked homozygous SNPs will be distributed around the phenotype causing SNP. Plotting the homozygous/heterozygous SNP ratio, should reveal a peak in the non-recombinant region, because of the high density of homozygous SNPs in that location (there is a consistent distribution of heterozygous SNPs from the mapping line across the genome). The causative SNP should be located at the peak of this ratio distribution.

### Creating a model genome

Running: ``ruby create_model_genome.rb dataset_name`` will generate a new model genome of that name based on Arabidopsis chromosome 4 and the experiment detailed above, in ``fragmented_genome_with_snps/arabidopsis_datasets/dataset_name``. This includes a FASTA file with the sequences of each fragment, and a VCF file with the SNPs on each fragment. In the INFO field of the VCF, each SNP has been given an allele frequency (AF). Heterozygous SNPs will generally have AF = ~0.5, and homozygous AF = ~1.0, but this will vary with pooled data. In the model, each SNP has been given an allele frequency of exactly 0.5 or 1.0. The ``dataset_name`` argument specifies the name of the directory to create and store the files in. In [create_model_genome.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/create_model_genome.rb), this currently places that directory: ``"~/fragmented_genome_with_snps/arabidopsis_datasets/dataset_name"`` (but this can of course be altered).

By editing [create_model_genome.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/create_model_genome.rb), the features of the model genome created by it can be altered. The variables ``hm_r`` and ``ht_r`` contain the R code needed to create the model homozygous and heterozygous SNP distributions respectively. The variable ``contig_size`` provides the minimum size for contigs, where the maximum size is double this value, and each contig's size is randomly chosen within this range. The ``fasta_file`` variable points to a FASTA format file with one chromosome, which in this case is *Arabidopsis thaliana* chromosome 4. This chromosome forms the basis of my model genome. Fragments (contigs) of the input FASTA file are created, and the SNPs are applied to them. A new FASTA file is created with the contigs randomly ordered, and a VCF file is created. These serve as the starting files for my algorithm. I have also included a number of other useful output files: the FASTA file with the contigs ordered correctly, txt files that list the homozygous and heterozygous SNPs, and another txt file (info.txt) that records the values of ``hm_r``, ``ht_r`` and ``contig_size``.

All the datasets I have created with [create_model_genome.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/create_model_genome.rb) are called ``10K_dataset...``,  see [arabidopsis_datasets](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/arabidopsis_datasets).

I have another script [small_test_genome.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/small_test_genome.rb) that can be used to make model genomes of length 2Kb, for quick testing of algorithm ideas. The datasets I have created with this, I have named ``small_dataset...``,  see [arabidopsis_datasets](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/arabidopsis_datasets).

Why not use real data straight away?
 - With model data, I know the location of the causative mutation and the correct order of the fragments. This means I can tell whether or not my algorithm works.
 - I can generate datasets that model different experiments/sample data quickly: e.g. different fragment sizes, different SNP distributions, different experiment

Re-ordering the Genome
----------

I aim to use a **genetic algorithm** to rearrange contigs, and find a permutation of the contig order that has the correct homozygous/heterozygous SNP ratio distribution. Permutations with distributions close to the expected, are likely be ordered in an arrangement close to the correct. The algorithm I am creating will eventually be capable of determining the location of the causative mutation, in a genome from the experiment detailed above. The algorithm will take a FASTA file with the unordered genome fragments (contigs), and a VCF file containing the SNP positions for each fragment.

### Genetic Algorithm

A genetic algorithm is a kind of iterative improvement algorithm. These are where an initial starting solution to a problem is incrementally improved, in order to reach an eventual perfect solution. In the case of a permutation problem, like my contig rearrangement, it is the *ordering* of elements (contigs) that needs to be changed until correct.

Genetic algorithms are based on the principles of natural selection. A "population" of solutions is created, and those with the highest "fitness" (correctness) are selected. New "offspring" solutions can be created through recombination of the components of "parent" solutions, and "mutations" can be introduced by creating novel components or parameters. If subsequent populations are created by the recombination of selected members of the current population, and/or the introduction of mutant solutions, the best solutions in each generation will be better than the best solutions from the previous generation. In this way, the perfect solution (in this case the perfect permutation of contig order) will eventually be reached.

With my model genome, I know the correct permutation of fragments/contigs (I have named each frag1, frag2...), so I am able to see if the algorithm has arranged them correctly (see 'Assesing the algorithms performance' below). I give each permutation of the fragment order created by the algorithm, a fitness score based on it's homozygous and heterozygous SNP distributions.

### Mutation

In my genetic algorithm, I use mutation methods to rearrange the contigs into novel permutations. The mutation methods used are from my [PMeth ruby gem](https://github.com/edwardchalstrey1/pmeth). The repo provides a description of what they do.
In each generation of my genetic algorithm, I create mutants of permutations from the previous generation.

### Permutation Fitness

Genetic algorithms require a fitness method, in order to tell how close each solution comes to solving the given problem. With a permutation problem, the fitness method should score the permutation based on how correct the ordering of elements is. I have tested several different fitness methods. Each of these makes a slightly different assumption about the SNP distributions of the genome.

For a detailed look at how the fitness methods work, see [Documentation/fitness](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/fitness_methods/fitness_methods.md)

### Quit condition

A genetic algorithm needs a way of knowing when to stop. Ideally, you want it to continue until a correct solution has been reached, but in reality this may be impossible. I have implemented a [quit condition](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Quit_condition/quit_condition.md) that terminates my genetic algorithm, when improvement (in fitness) is judged to be negligible.

How my genetic algorithm works and how to run it
----------

The ``evolve`` method in the [GATOC class](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) is a wrapper method that forms the bulk of my genetic algorithm. A script is required to call ``GATOC.evolve``, and whilst testing the algorithm I have been using [algorithm_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_test.rb).

To run [algorithm_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_test.rb): ``ruby algorithm_test.rb dataset_name run``, where ``dataset_name`` is the same as the one used to create the model genome, and ``run`` is the name of the directory in which to save results. Subsequent required command line parameters are described in [algorithm_test.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/algorithm_test.rb): these correspond to input parameters for ``GATOC.evolve``.

When running the algorithm with real data, a new ruby script will be needed to call ``GATOC.evolve``

For a detailed look at how the algorithm works see [Documentation/evolve](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/evolve.md)

Assesing the Genetic Algorithm's performance
---------------------------------------

My experiment relies on the premise that the fitness of a permutation, as I have calculated it, is related to how close that permutation is to the correct arrangement of contigs. To test this, I use a number of distance metrics, which compare the order of permutations and provide a similarity score. I have made a ruby gem called [PDist](https://github.com/edwardchalstrey1/pdist) that contains the methods for calculating these scores. Descriptions of how they work can be found in the repository.

### Testing the algorithms fitness methods, and the distance metrics from PDist

I have carried out an experiment to determine which of the fitness methods I have come up with are most worth pursuing, and which distance metrics (from  [PDist](https://github.com/edwardchalstrey1/pdist)) are most useful for analysing my results. For my findings, see [Progress/Do_metrics_work](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Do_metrics_work/do_metrics_work.md).

### Data analysis

For information on how I analysed the output from my algorithm see [Documentation/analysis](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/analysis/analysis.md)

### Candidates for the causative SNP mutation

For information on how I worked out candidates for the causative SNP mutation, see [Documentation/causal_mutation](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Documentation/causal_mutation/causal_mutation.md)

Results
-----

Below are the results of running many replicates of my algorithm with different fitness methods.

[Results for count_ratio fitness method](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Results1_count_ratio/results.md)

[Results for snp_distance fitness method](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Results2_snp_distance/results.md)

[Results for max_density fitness method](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Results3_max_density/results.md)

[Results for max_ratio fitness method](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Results4_max_ratio/results.md)

[Results for max_hyp fitness method](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/Progress/Results5_max_hyp/results.md)

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
