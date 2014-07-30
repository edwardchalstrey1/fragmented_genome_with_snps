Locating the causal mutation, whilst testing my algorithm on my model genome
========================================================

Script: [find_causal_mutation.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/find_causal_mutation.rb)

Running ``ruby find_causal_mutation.rb dataset run last_generation`` will print out the position of the homozygous SNP, at the centre of the peak in homozygous to heterozygous SNP ratio, for both the ordered genome and the fittest permutation outputted from the genetic algorithm. ``dataset`` is the directory within ``fragmented_genome_with_snps/arabidopsis_datasets`` that contains the input data (FASTA and VCF files) for the algorithm. ``run`` is the directory within ``dataset`` that contains the output from one run of the algorithm. ``last_generation`` should be an integer of the final generation in the run.

Locating the causal mutation from real world data (and how the location works)
============

After running my genetic algorithm with real world data, the causal mutation could be found using a modified version of ``find_causal_mutation.rb``. Use the methods ``find_peak`` and ``closest_snp`` drom the [LocateMutation](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/locate_mutation.rb) class as shown below:

```ruby
peak =  LocateMutation.find_peak(hyp, n) # Find the peak in the approximated (hypothetical SNP) distribution
causal_mutation = LocateMutation.closest_snp(peak, hm)
```

``hyp`` is a hypothetical distribution that represents the ratio of homozygous to heterozygous SNPs (see [SNPDist.hyp_snps](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/lib/snp_dist.rb)) across the genome/permutation.
``n`` is the number of equally spaced points at which the density of ``hyp`` SNPs are estimated, and should be a power of two.
``hm`` is the distribution of homozgous SNPs across the genome/permutation, that was used with the heterozygous distribution to create ``hyp``.

The ``causal_mutation`` therefore, is the homozgous SNP that lies closest to the highest peak in the ratio of homozgous to heterozgous SNP density.