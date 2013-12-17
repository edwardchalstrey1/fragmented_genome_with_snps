Rearrangement method ideas
========================================================

Machine learning
-----

Learning involves giving the machine examples, then the more examples, the better the machine gets at telling how a new thing should be catagorized. For example, a machine learning system could be trained on email messages to learn to distinguish between spam and non-spam messages.

Hill Climbing algorithms
------

An arbitrary solution to a problem is generated, then an element of the solution is changed, if this produces a better solution, it is kept. Iteratively, the solution gets better. It can however, reach a local maxima, where changing one element cannot make a better solution, whilst the global maxima (the best possible solution) is not actually reached.

What to do
--------

Looking at the [NGM website](http://bar.utoronto.ca/ngm/cgi-bin/emap.cgi), with sample data, we can see that there is a noticeable homo/hetero signal distribution around the causative mutation.

You need to modify your model genome:

- Two sets of SNPs (hetero and homozygous), that will generate a normally distributed (or similar curve to NGM sample) ratio distribution
- The ratio is maximized at the causative mutation, and distributed around it
- The standard deviation/ curve of the distribution is down to the kernel density. How tightly does this fit around the data? You decide: we want the distribution to look as normal as possible, as this will be what we try and regenerate when putting the fragments back together
- Use a qq plot to make sure the ratio is a normal distribution
- then use the qq plot to test the distribution of rearranged fragments' hom/het ratio (score the rearrangements based on how close to the norm dist they are

Questions
========================================================

Back-crossing experiment produces F2 EMS mutants, with SNPs distibuted around causative mutation.
Sequencing and contig assembly produces fragments with SNPs at known positions.

1. What will the exact function of my tool be?
 - Locate the causative mutation by looking at the fragments
 - The SNPs will form a homozygous/heterozygous ratio that is normally distributed about the causative mutation?

2. How should I model this? 
 - If I put "homozygous" and "heterozygous" SNPs into the fasta/vcf data, how will I then re-identify them? NGM uses chastity statistic to measure the frequencies of hom and het SNPs from pooled data. But the experiment I'm modelling is different? If not, then will I not just be replicating NGM?
 - I have come up with a way of making a hom/het ratio distribution, [see ratio.R](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/ratio/ratio.R)
 - Once I have made new datasets for this version of the model genome, I can read off generated fasta/vcf files and put frags in correct order and write some code that can work out/display the ratio distribution. This code will be used later when rearranging randomly ordered fragments.
 - It might be better to just use SNP distribution, rather than ratio distribution: however, as we have seen, the SNPs in the recombinant region aren't neccesarily normally distributed, this is why we went for ratio in the first place

3. How will the algorithm work?
 - I want to read off the data e.g. fasta/vcf and regenerate the normal distributed ratio, which will be done by putting the fragments back into the correct order
 - The fragments will be rearranged by a hill climbing algorithm: this will iteratively attempt to reorder the fragments in a way that creates a hom/het SNP ratio normal distribution (fragments with zero ratio need not be ordered correctly/ may be ignored completely)
 - The causative mutation will be at the centre of the distribution: see figure source [see ratio.R](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/ratio/ratio.R). The distribution is reformed by reading the SNP positions from the fragments, and working out their SNP position in the genome (for this fragment order). The SNP lists, hom and het, can then be used to work out the ratio distribution.


Hill Climbing
---------

1. **A reasonable starting solution is required:**
For my fragment rearrangement algorithm I propose using my skew/gradient method. Even if this method proves to be of equal uselessness for reforming the narrow ratio distribution, than the wider SNP distribution I have tested it on so far, It can at least split the fragments into 2 groups, and then they can be ordered by SNP density. One problem I may come across is that ordering the fragments by SNPs per nucleotide (density) won't be very useful for low SNP numbers (within the "non-recombinant region" where there are just homozygous SNPs). However, the low number of SNPs does mean that the normal distribution will almost definitely only be reformed by a correct rearrangement.

2. **Evaluate correctness funtion:**
This would be an extension of the code I came up with earlier that calculates the ratio distribution (for the fragment rearrangement), but would score it on how normal it is. I could use a qq plot to do this.

3. **Operator to modify complete solution:**
NOT SURE YET.
Thinking about it, for this to work, we will want to deal with all of the fragments, some of them from outside the recombinant region may have a noticeable low sub-peak of ratio (i.e. bumps in the overall ratio plot). To make the model genome more realistic, you may want to add some more of these "bumps".

- Hopefully the initial solution will reveal some distinct peaks

Most appropriate type of hill climbing algorithm
-----

Hill climbing may be limited for fragment rearrangement if all I do is change one fragment's position at a time. There will be lots of local maxima rearrangement states.

### Genetic search algorithm/ evolutionary computation - type of beam search
Could use multiple starting states, e.g. the current rearrangement methods, or just random ones. Their normality is determined with qq plot (or whatever). Then the most normal ones get to produce offspring rearrangements. "Genes" that get passed on and recombined can be X fragments, the order within each X conserved, but each time we do a "cross", the value of X should be different, so wrongly ordered fragments don't get stuck together. We can also introduce random "mutations", one of the "genes" of X fragments could be shuffled/reversed. The "fittest" (most normal) rearrangements should be crossed in each generation, but this should include the parents, in case all of the children are worse.
