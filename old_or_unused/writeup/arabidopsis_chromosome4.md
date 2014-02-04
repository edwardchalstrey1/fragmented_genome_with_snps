*Arabidopsis thaliana* chromosome 4, with SNPs
========================================================

A more realistic model genome
--------

To better evaluate the effectiveness of the skew method (method 2b) of fragment rearrangement (see [rearrangement methods](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md)), I have created a new model genome ([arabidopsis_c4_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_c4_w_snps.rb)), from *Arabidopsis thaliana* chromosome 4 data available [here](ftp://ftp.arabidopsis.org/Sequences/whole_chromosomes/) (TAIR10_chr4.fas). I have given the chromosome 70,000 SNPs normally distributed about a causative mutation (so 70,001 total), which is located at position 10,000,000 of the 18,585,056 base chromosome (golden path length). The standard deviation of the distribution is 2,000,000. 70,000 was the chosen number of SNPs: searching online and finding [this paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3290786/) which describes the number of SNPs between two ecotypes of *Arabidopsis* as around 350,000, and dividing by 5 (the number of chromosomes). The exact number of SNPs doesn't matter for this model genome, but this working was just to get the order of magnitude appropriate.

The fragment (contig) sizes are 10-20Kb.

I have updated [arabidopsis_c4_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_c4_w_snps.rb) and [json->fasta.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/json-%3Efasta.rb), so that new datasets (new fasta and vcf files) for the arabidopsis chromosome with SNPs can be generated relatively quickly, with different SNP positions (and fragments) each time, but the same number of SNPs and the same causative mutation.
For this model genome I cannot perform the shapiro wilk test normailty because the test uses a maximum 5K sample size. This no longer matters because the SNP positions are taken directly from the rnorm function in R, and therefore definitely follow a normal distribution. See [arabidopsis_c4_w_snps.rb](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_c4_w_snps.rb) normal_dist method for details.

It is difficult to say at this stage how *realistic* this version of the model genome is, but it should provide me with a better basis for testing methods of fragment re-arrangement, and eventually causative mutation location.

For more information, and what I mean by *realistic* see the repository [README](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/README.md)




