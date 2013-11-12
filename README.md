Fragmented Genome with SNPs
===========================

I am creating a model genome (or chromosome/sequence) from an individual created by the out-crossing of a mutant individual with non-mutants, where the mutation is an experimentally induced, phenotype altering SNP. Specifically what this models is the following out-crossing experiment: a beneficially mutated individual is created by experimentally inducing SNPs (e.g. using a chemical mutagen on a plant to create a disease resistance mutation). The phenotype altering SNP needs to be identified. The mutant individual is out-crossed (bred) with a non-mutant individual that has no experimentally induced SNPs. Out-crossing like this multiple times, but always selecting for progeny with the mutant phenotype, will result in an individual that has a region of high SNP density near the mutant position, due to linkage, and a much lower number of SNPs in other parts of its genome. The remaining SNPs are normally distributed around the causative mutation, because with each out-cross, the SNPs closest to mutation being selected for have a higher chance of being conserved through linkage, than SNPs further from the location of the mutation.

After creating a model genome based on the expected result of the out-crossing experiment detailed above, I am designing an algorithm that will locate the causative mutation, based on the distribution of SNPs. To emulate real data, I first split the genome into fragments, which model contigs assembled from high-throughput sequencing reads. My algorithm will need to rearrange these fragments into their proper order, then locate the causative mutation.

For a chronological writeup of my results of algorithm development (and model genome improvements), see the files in the writeup folder in this order: [Rearrangement methods](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md), [model genome dataset 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/dataset2.md), [rearrangement methods part 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/p2_rearrangement_methods.md), [model genome dataset 3 and subsequent datasets](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/arabidopsis_chromosome4.md), rearrangement methods part 3...

**Details of the main files in the repo are below. For a more comprehensive summary of how/why they are used, read the writeup**

### Modeling a fragmented genome with SNP's in Ruby, creation of dataset 1 and 2:

The mutant position is located at position 100,000 in [dataset 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/fasta_vcf) and [2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/fasta_vcf_d2)), the mid-point of a 200Kb sequence. The sequence is split into fragments of 50-250b length.

1. [fragments_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/fragments_w_snps.rb) - This generates the fragments described above. JSON files are created: one contains the fragments, another the SNP positions for each fragment.

2. ^ a text file of the snp positions is also created. Running [normality_test.txt](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/normality_test.txt) in R confirms that the snp positions are normally distributed:
 - see [dataset2.md](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/dataset2.md) for details

3. [json->fasta.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/json-%3Efasta.rb) - uses the information in the JSON files to construct a rudimentary fasta format file and VCF file. (see commits from before 11/11/13 for dataset 1 and 2)

4. [dataset2.md](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/dataset2.md) contains details of the improvements made to the model genome after testing [fragment rearrangement methods](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/rearrangement_methods.rb), see below

### A more realistic model: *Arabidopsis thaliana* chromosome 4

1. Dataset 3 and all subsequent datasets are based on [*Arabidopsis thaliana* chromosome 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/arabidopsis_chromosome4.md).
 - [arabidopsis_c4_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_c4_w_snps.rb) replaces [fragments_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/fragments_w_snps.rb)
 
2. When running [arabidopsis_c4_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_c4_w_snps.rb) and [json->fasta.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/json-%3Efasta.rb), do as follows for a named dataset e.g. datasetX
 - ruby "" datasetX
 - the name of the dataset entered into the command line saves the outputted files in an appropriate location

### Designing an algorithm that will determine the position of a phenotype altering mutation:

1. [density_method_testa.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/density_method_testa.rb) - uses information from fasta and VCF files to work out the SNP density of each fragment in the fasta file, measured as SNPS per Kb.

2. Rearranging the fragments into their original order:
 - [rearrangement_methods.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/rearrangement_methods.rb) contains rearrangement methods
 - [see rearrangement_methods.md]() for details of methods
 - then see [dataset_2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/dataset2.md) and [part 2 of rearrangement methods](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/p2_rearrangement_methods.md)
 
3. Identifying the causative mutation
