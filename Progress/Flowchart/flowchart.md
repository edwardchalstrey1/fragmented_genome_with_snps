Genetic algorithm flowchart
========================================================

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/evolve.png?raw=true)

### FASTA file to population

**(1-2)** [ReformRatio::fasta_array](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) obtains A FASTA object for each contig, and adds them to an array. [GATOC::initial_population](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) randomly orders the FASTA object array, and adds it to a super-array (the population). This is repeated to get the population of the desired size.

### VCF file to SNP positions on each contig

**(3-4)** [ReformRatio::get_snps](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) creates arrays of the SNP positions on each fragment, and outputs arrays and hashes of the other information in the VCF. *All this information is stored in memory.*

### SNP positions for a contig permutation

**(5-6)** [ReformRatio](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) contains other methods that obtain arrays of the homozygous and heterozygous SNPs across the entire genome, for each contig permutation. *Perhaps these methods can be grouped into one method, to avoid large arrays of the data being stored in memory*

### Calculating the ratio of homozygous to heterozygous SNPs, for each permutation

**(7-8)** [FitnessScore::count](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/fitness_score.rb) is used to get the number of SNPs at a set number of divisions across the genome. This is done for homozygous and heterozygous SNPs. Then [FitnessScore::ratio](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/fitness_score.rb) creates an array of the ratio of homozygous to heterozygous SNPs for each of the genome divisions.

### Permutation fitness

**(9-10)** [FitnessScore::score](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/fitness_score.rb) calculates the "fitness" of a permutation, based on the comparison of it's distribution of homozygous to heterozygous SNP ratio, with an expected distribution. Currently (June 25th) I am using the distribution from correctly ordered contigs as the expected distribution, since I am still testing whether the genetic algorithm aproach will work.
THEN SAY HOW METHOD WORKS
