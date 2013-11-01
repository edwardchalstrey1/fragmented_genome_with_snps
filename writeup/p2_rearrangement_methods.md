Rearrangement methods part 2
========================================================

Having created a new dataset ([dataset2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/dataset2.md)), I will re-run the existing fragment rearrangement methods described in [rearrangement_methods.md](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md). I will then continue to write up additional experiments here. I used this [fasta](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/fasta_vcf_d2/frags_shuffled.fasta) and this [VCF](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/fasta_vcf_d2/snps.vcf). Whilst dataset 1 contained 1310 fragments, dataset 2 contains 1321. 

Existing Method Scores
----------------------

Below are a list of the methods used in [rearrangement_methods.md](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md), and a table comparing the scores for dataset 1 and 2:

Highest possible score (C0), Density order Score (C1), Random Score (C2), Even Odd Method Score (M1a), Odd Even Method Score (M1b), Left Right Method Score (M2a), Left Right Density Method Score (M2b)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/dataset_scores_table.png?raw=true)
[Figure 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/dataset_scores_table.png)

As you can see from Figure 1, for each rearrangement method, the ordinal similarity score is similar in both datasets (and slightly higher for dataset 2 in most cases, probably due to the larger number of fragments). This was what I expected for the controls and method 1 (odd/even), since the change in the model genome was aimed at maximising the potential of method 2b (left right density), see [dataset2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/dataset2.md). This method relies on the skew of SNPs to determine each fragments position in the rearranged order, see [rearrangement_methods.md](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md).

Figure 1 also shows that method 2b has a similar ordinal similarity score for both datasets. What this may suggest is that the skew of SNPs on a fragment is in fact not a good indicator of its position in the genome. However, it could be that the current skew method (2b) is simply incorrect in its placing of the fragments.

For SNP density plots from dataset 2 (and bar plot of scores) see [dataset 2 plots.](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/figures/figures_d2)

To investigate how each fragments' skew relates to it's position in the proper order, I have plotted graphs of the "SNP gradient" for each fragment. The R script used to create figure 2, 3 and 4 can be found [here](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/skew_scatta.R). 
My prediction was that the gradients of frags would show a similar pattern to the normal distribution of SNPs, with the gradients being higher towards the midpoint of the model genome (where the mutation is), and lower at either end. However, I also thought that some of the fragments right in the centre of the genome might have a lower gradient, due to their having many SNPs across their lengths (i.e. on either side, reducing the skew). In figure 2 below, I use an absolute value of gradient. It will also be useful to show the difference between the positive and negative gradient values,  to determine whether the gradient (and skew) can be used to tell the side of the model genome (either side of the mutant point) that a fragment is from (figure 4).

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/skew_scatter.png?raw=true)
[Figure 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/skew_scatter.png)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/skew_scatter2.png?raw=true)
[Figure 3](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/skew_scatter2.png)
![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/example_gradient_f629.png?raw=true)
[Figure 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/example_gradient_f629.png)

Figures 2 and 3 are plots of each fragment in their proper order, against their gradients. In figure 2 the gradients used are the absolute values of the gradients used in figure 3. 

The "SNP gradients" are determined by creating plots of the nucleotides against the number of SNPs for each fragment, where the number of SNPs is a value of 0 or 1.
Figure 4 shows an example plot, where the gradient of the blue line (a linear model fitted to the data), is the gradient for that fragment plotted in figure 2/3. For each fragment, the skew of SNPs to either side of the fragment produces a positive gradient (right skew) or negative gradient (left skew, as with fragment629 in figure 2). Calculating the skew in a way that could be quantified was important, in order to plot figure 2/3.

Interpreting what figure 2 and 3 show about the SNP gradient/skew of fragments is difficult, because of the large number of fragments that have a gradient of zero (or close to zero). As you can see, the loess smoothing lines (red) for figure 2 and 3 are drawn across value zero of gradient. The reason for this is that many of the fragments have either 0 SNPs (which gives them a gradient of exactly zero) or 1 SNP, or another low number of SNPs. 

This is another problem with my model genome that I have identified; the gradients and skew may reveal more about the position of fragments if a larger number of SNPs is used. Fragments with just 1 SNP (or any low number) are less likely to have a useful skew/gradient than fragments with a larger number of SNPs. If the SNPs in my model genome are not very concentrated with regard to the fragment size, the position of SNPs on individual fragments could be a lot more arbitrary. If the fragments were larger (in proportion to the genome size) and/or the genome had a higher number of SNPs around the mutant position, the fragments may exhibit a much clearer skew, that is related to their position in the genome. Figure 5 below illustrates this point.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/example_gradient_f321.png?raw=true)
[Figure 5](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/example_gradient_f321.png)

Working out the gradient in the way I have done may not be useful for fragment rearrangement. However, it does show that the skew method (method 2b) is not likely to be effective with the model genome used to create datasets 1 and 2.

Dataset 1 and 2 are from a model genome of 200Kb in length with 1000 SNPs normally distributed around the mutant position 100,000, and fragment sizes of 50-250b.

