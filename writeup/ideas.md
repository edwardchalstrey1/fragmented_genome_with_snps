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

Info
========================================================

Back-crossing experiment produces F2 EMS mutants, with SNPs distibuted around causative mutation.
Sequencing and contig assembly produces fragments with SNPs at known positions.

1. What will the exact function of my tool be?
 - Locate the causative mutation by looking at the fragments
 - The SNPs will form a homozygous/heterozygous ratio that is normally distributed about the causative mutation

### How should I model this? 
1. If I put "homozygous" and "heterozygous" SNPs into the fasta/vcf data, how will I then re-identify them? NGM uses chastity statistic to measure the frequencies of hom and het SNPs from pooled data. But the experiment I'm modelling is different? If not, then will I not just be replicating NGM? No:
 - I have come up with a way of making a hom/het ratio distribution, [see ratio.R](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/ratio/ratio.R)
 - I can use the VCF ruby module to assign hom and het SNPs different frequency values or something
 - Once I have made new datasets for this version of the model genome, I can read off generated fasta/vcf files and put frags in correct order and write some code that can work out/display the ratio distribution. This code will be used later when rearranging randomly ordered fragments.
- It might be better to just use SNP distribution, rather than ratio distribution: however, as we have seen, the SNPs in the recombinant region aren't neccesarily normally distributed, this is why we went for ratio in the first place
- Multiple datasets: different SNP ditributions (use the one in ratio.R first), and then also different fragment sizes (less important)
- Try and make the SNPs as realistic as possible
- Allele frequency = what proportion of the sample has the ALT. So for just one model genome, heterozygous SNPs AF=0.5 and homozygous AF=1.0

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

Using Brute Force
============

1. If there are X fragments, how many rearrangements are there? 
 - There will be X! permutations. See [permutations without repetition](http://www.mathsisfun.com/combinatorics/combinations-permutations.html). ! means the factorial function.

2. If the number of fragments, and by extension X!, is low enough, can we compute the correct rearrangement in acceptable time?
 - If so, this would eliminate the need for hill climbing or genetic algorithm
 
Plan of action
===============

QQ plot
-------

I will use a QQ plot to score the SNP ratio distribution against a random normal distribution. This will form the basis by which rearrangements are ranked, by giving each a score value.

Implementing a genetic algorithm in ruby
--------------------

[Potential drawback](http://www.solver.com/genetic-evolutionary-introduction) of genetic algorithm

**See THE CROSSOVER page of lab book for the crossing/recombination method**

Genetic algorithm in R
---------

### [GA package R](http://cran.r-project.org/web/packages/GA/GA.pdf) - [GA paper](http://www.jstatsoft.org/v53/i04/paper)

**Definitions**

1. Search space: all the possible solutions, each soultion is one point in the search space

**Aguments**

1. type: "permutation" (involves re-ordering a list)
2. fitness: an R function, of fitness for the rearrangement/permutation, QQ plot. The function must have a string input, which represents the solution (rearrangement), and returns a number (fitness score)
 - ... : additional arguments for fitness function, allows some variables for the fitness function to be fixed during the search
3. min: a vector of the decision variables that give the minimum of the search space **NOT ENTIRELY SURE WHAT THIS MEANS**
4. max:
5. nBits: *not needed*
6. population: R function for randomly generating a set of fragment rearrangements (initial population)
7. selection: R function that creates a sub-population of the population, based on a threshold of the fitness (i.e. selecting those above a normality threshold)
8. crossover: R function that crosses rearrangements. It should split the permutation into chunks of X size
9. mutation
10. popSize
11. pcrossover
12. pmutation
13. elitism
14. monitor
15. maxiter
16. run
17. maxfitness
18. names
19. suggestions
20. keepbest
21. parallel
22. seed
