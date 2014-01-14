Fragmented Genome with SNPs
===========================

I am creating a model genome (or chromosome/sequence) from an individual created by the out-crossing of a mutant individual with non-mutants, where the mutation is an experimentally induced, phenotype altering SNP. Specifically what this models is the following out-crossing experiment: a beneficially mutated Arabidopsis individual is created by experimentally inducing SNPs with EMS (e.g. a disease resistance mutation). The phenotype altering SNP needs to be identified. The mutant individual is out-crossed (bred) with a mapping line that has no experimentally induced SNPs. This mapping line is a cross between two Arabidopsis ecotypes (one of which the same ecotype as the mutant), and therefore has heterozygous SNPs across its genome. Out-crossing like this multiple times, but always selecting for progeny with the mutant phenotype, will result in an individual that has a non-recombinant region with high homozygous to heterozygous SNP ratio near the mutant position, due to linkage. The homozygous SNPs are distributed around the causative mutation, and with each out-cross, the SNPs closest to mutation being selected for have a higher chance of being conserved through linkage, than SNPs further from the location of the mutation. Plotting the homozygous/heterozygous SNP ratio should reveal the non-recombinant region with a normal distribution of homozygous/heterozygous SNP ratio, the causative SNP should be located at the peak of this ratio distribution.

After creating a model genome based on the expected result of the out-crossing experiment detailed above, I am designing an algorithm that will locate the causative mutation, based on the distribution of SNPs. To emulate real data, I first split the genome into fragments, which model contigs assembled from high-throughput sequencing reads. My algorithm will need to rearrange these fragments into their proper order, then locate the causative mutation.

Creating a model genome
--------------------

Running: **ruby [arabidopsis_c4_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_c4_w_snps.rb) dataset_name** will generate a new model genome of that name based on Arabidopsis chromosome 4 and the experiment detailed above, in fragmented_genome_with_snps/arabidopsis_datasets. This includes a FASTA file with the sequences of each fragment, and a VCF file with the SNPs on each fragment. In the INFO field of the VCF, each SNP has been given an allele frequency (AF). Heterozygous SNPs will generally have AF = ~0.5, and homozygous AF = ~1.0, but this will vary with pooled data. In the model, each SNP has been given an allele frequency of exactly 0.5 or 1.0.

Locating the causative mutation
--------------

I am currently attempting to use a genetic algorithm to rearrange the fragments, by reforming the homozygous/heterozygous SNP ratio distribution.

1. [reform_ratio.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/ratio/reform_ratio.rb) uses my own implementation of a genetic algorithm. To run it: **ruby reform_ratio.rb dataset_name**, where the dataset name is the same as the one used to create the model genome.

2. [charlie_1.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/ratio/charlie_1.rb) is an attempt to use the [charlie](http://charlie.rubyforge.org/) ruby gem genetic algorithm. To run it: **ruby charlie_1.rb dataset_name**.
 - currently, the algorithm is not rearranging the frags that the fitness function is acting on, so the fitness value being displayed is random

3. [ratio.R](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/ratio/ratio.R) is used within the fitness function of my genetic algorithm. It compares the homozygous/heterozygous SNP distribution of the rearranged fragments, to the same distributions I used when creating the model genome, using a qq plot. A coefficient value is obtained, between 0 and 1. The closer the value is to 1, the more likely it is that the correct fragment arrangement has been found.
