Genetic algorithm flowchart
========================================================

<script src="http://www.gliffy.com/diagramEmbed.js" type="text/javascript"> </script><script type="text/javascript"> gliffy_did = "5899757"; embedGliffy(); </script>

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/genetic_algorithm.png?raw=true)

### FASTA file to population

**(1-2)** [ReformRatio::fasta_array](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) obtains A FASTA object for each contig, and adds them to an array. [GATOC::initial_population](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/GATOC.rb) randomly orders the FASTA object array, and adds it to a super-array (the population). This is repeated to get the population of the desired size.

### VCF file to SNP positions on each contig

**(3-4)** [ReformRatio::get_snps](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) creates arrays of the SNP positions on each fragment, and outputs arrays and hashes of the other information in the VCF. *All this information is stored in memory.*

### SNP positions for a contig permutation

**(5-6)** The [ReformRatio](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/reform_ratio.rb) class contains other methods that obtain arrays of the homozygous and heterozygous SNPs across the entire genome, for each contig permutation. *Perhaps these methods can be grouped into one method, to avoid large arrays of the data being stored in memory*
