Rearrangement methods part 3
========================================================

Having updated the genome model being used to [*Arabidopsis thaliana* chromosome 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/arabidopsis_chromosome4.md), I will continue to evaluate my methods of fragment rearrangement here.

To re-cap: I need create an algorithm that re-orders the fragments of the genome into their proper order, then I will go on to include this algorithm in my algorithm for locating the causative SNP mutation. See the [README](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/README.md) for details.

Current rearrangement methods
-----------

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/scores_table.png?raw=true)
[Table 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/scores_table.png)

The method ids in figure 1 are as follows: Highest possible score (C0), Density order Score (C1), Random Score (C2), Even Odd Method Score (M1a), Odd Even Method Score (M1b), Left Right Method Score (M2a), Left Right Density Method Score (M2b). The methods themselves are described in [rearrangement methods part 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md).

I have generated 4 separate datasets (model genomes) based on the arabidopsis chromosome, each with different SNP positions and fragments, but based on the same random normal distribution of positions, and same fragment size range. See [*Arabidopsis thaliana* chromosome 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/arabidopsis_chromosome4.md) for details. In table 1, the means of the scores for these datasets are shown in the 'datasets 3-6' column. Whilst dataset 1's genome was split into 1310 fragments and dataset 2 1321 fragments, the average number of fragments for these four arabidopsis datasets is approximately 1236. This explains the lower scores for datasets 3-6, compared with 1 and 2. For details on how the score is calculated, see [rearrangement methods part 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md). The lower the score, the closer the rearranged fragment order is to the proper order.

As you can see from table 1, the existing rearrangement methods have high scores for datasets 3-6, and improve little upon the controls. This was what I expected for most of the methods, however I did expect method 2b to perform better, due to the higher number of SNPs per fragment for these datasets. Method 2b: left right density method [(read decription of method 2a/b in rearrangement_methods.md)](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md) uses the fragment's SNP skew, and their SNP densities (measured as SNPs per Kb) to rearrange them. As I explained in [rearrangement methods part 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/p2_rearrangement_methods.md), I need to evaluate how well a fragment's SNP skew, and density, can be used to determine it's position in the model genome.

SNP density as a single value for each fragment, measured as SNPs/Kb
------------------

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/d_o.png?raw=true)
[Figure 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/figures/figures_d2/d_o.png) shows SNP density for each fragment in **dataset2** 

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/d_o.png?raw=true)
[Figure 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/d_o.png) shows SNP density for each fragment in **dataset6**

Figures 1 and 2 show that the SNP density of fragments in the arabidopsis model genome (dataset6) are significantly closer to a true normal distribution, than those in the previous model genome (dataset2). This is because the older model genome had a large proportion of fragments with zero, or close to zero SNPs.

The reason that the SNP density in figure 2 does not show a perfect normal curve, despite the SNP positions in my model dataset being generated with the rnorm function in R, is that the fragments are not of equal length. Since this is also the case with real contigs, it follows that a real dataset would show a similar pattern to figure 2, if we assume that the SNP positions in the real data are also normally distributed. Since this is what I *am* assuming, I think that the model genome based on arabidopsis chromosome 4 that I have created is a fairly realistic version of what I have described in the project [README](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/README.md).

Figure 1 and 2 were generated with this [ruby script](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/scatter_graphs.rb), which uses the function created in this [R script](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/scatter_vectors.R). See [here](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/tree/master/arabidopsis_datasets/dataset6/figures) for all the figures it created for dataset 6, one for each rearrangement method, the same as was done for dataset 1 in [rearrangement methods part 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md). The figures for each method are marked by their method ids (see figure 1) in the filename.

Testing the skew for functionality
----------

It may be that trying to rearrange all the fragments based on SNP skew is ineffective, because the skew is only useful in determining the position of fragments with a high enough SNP density: where the skew determines the fragments place in the model genome, based on the normal distribution of SNPs, [see rearrangement methods](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md) - method 2 for details.

To evaluate how the skew of SNPs within fragments is related to their position in the genome, I have again worked out the gradients of SNP density within each fragment. Instead of the method I used for gradient determination in rearrangement methods part 2, I have used a new method that works out the kernel density estimate.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad.png?raw=true)
[Figure 3](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs.png?raw=true)
[Figure 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f257.png?raw=true)
[Figure 5](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f257.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f687.png?raw=true)
[Figure 6](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f687.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f1042.png?raw=true)
[Figure 7](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f1042.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/genome_kdens_snps.png?raw=true)
[Figure 8](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/genome_kdens_snps.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad_thresh.png?raw=true)
[Figure 9](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad_thresh.png)

What figures 9 an 10 suggest is that fragments further from the centre of the SNP distribution are likely to have a slightly higher skew (steeper gradient) than those closer to the centre, but that the gradient can be positive or negative.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs_thresh.png?raw=true)
[Figure 10](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs_thresh.png)

