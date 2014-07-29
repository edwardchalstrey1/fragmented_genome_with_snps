Locating the causal mutation, whilst testing my algorithm on my model genome
========================================================

Running ``ruby find_causal_mutation.rb dataset run last_generation`` will print out the position of the homozygous SNP, at the centre of the peak in homozygous to heterozygous SNP ratio, for both the ordered genome and the fittest permutation outputted from the genetic algorithm. ``dataset`` is the directory within ``fragmented_genome_with_snps/arabidopsis_datasets`` that contains the input data (FASTA and VCF files) for the algorithm. ``run`` is the directory within ``dataset`` that contains the output from one run of the algorithm. ``last_generation`` should be an integer of the final generation in the run.
