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

To evaluate how the skew of SNPs within fragments is related to their position in the genome, I have again worked out the gradients of SNP density within each fragment. Instead of the method I used for gradient determination in [rearrangement methods part 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/p2_rearrangement_methods.md), I have used a new method that works out the kernel density estimate. 

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad.png?raw=true)
[Figure 3](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs.png?raw=true)
[Figure 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs.png)

As you can see from figure 3 and 4, neither the gradient or absolute gradient are showing much of a pattern when plotted against the fragments in their original order. I expected there would be a much larger difference in the gradients of fragments closer to, and further from, the centre of the normal distribution of SNPs. Understanding this difference could have eliminated certain fragments, when designing the algorithm that searches for the causative mutation.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f257.png?raw=true)
[Figure 5](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f257.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f687.png?raw=true)
[Figure 6](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f687.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f1042.png?raw=true)
[Figure 7](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/example_gradient_f1042.png)

Figures 5-7 are examples of the gradient determination step used to generate figures 3 and 4. One inportant consideration here is that I have only included fragments with 6 or more SNPs. Fragments with a low number of SNPs are less likely to have a useful skew/gradient than fragments with a larger number of SNPs. This is because the position of SNPs could be a lot more arbitrary, and less as a result of the SNP distribution across the genome. The kernel density estimate is plotted for the SNPs across each fragment, and then the coef function is used to extract the coefficient of the linear model of x and y. See the R pseusocode below for details.


```r
# snp_pos <- vector of SNP positions for a fragment kernel_density_est <-
# density(snp_pos, n=((fragment_length)/5)) x <- kernel_density_est$x y <-
# kernel_density_est$y gradient <- coef(lm(y ~ x))
```


In the density function, n is the number of equally spaced points at which the density is to be estimated. To make the size of the sections that this will split each fragment into consistent, I divide the fragment length by a constant, 5. This means that for each fragment, the density estimate will be determined approximately for each 5 bases. 

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/genome_kdens_snps.png?raw=true)
[Figure 8](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/genome_kdens_snps.png)

Since I have used the kernel density estimate to determine the SNP gradients, I have plotted the kernel density estimate for the entire model genome (figure 8). What this shows is that the SNP density, as worked out by the kernel density estimate, is an accurate way of viewing the normal distribution of SNPs. Figure 8 also leads me to believe that the fragments that cover the steeper parts of the distribution, should show a skew of SNPs toward the centre (of the distibution).

Looking back at figures 3 and 4, I noticed that most of the fragments had a gradient of below 5.0e-09 (absolute), so I decided to plot these graphs again, excluding the fragments with gradients above that threshold.

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad_thresh.png?raw=true)
[Figure 9](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_grad_thresh.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs_thresh.png?raw=true)
[Figure 10](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs_thresh.png)

What figures 9 and 10 may suggest, is that fragments further from the centre of the SNP distribution are likely to have a slightly higher skew (steeper gradient) than those closer to the centre, but that the gradient can be positive or negative. However, eliminating the outliers like this shouldn't be neccesary to observe a true pattern.

To investigate further if the skew is useful for fragments with a high number of SNPs, I have eliminated fragments with low SNP numbers:

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs30.png?raw=true)
[Figure 11](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs30.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs50.png?raw=true)
[Figure 12](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset6/figures/skew_scatter_abs50.png)

I think that figure 11 and 12 illustrate clearly that the fragment's SNP skews will not be useful in rearranging them to the correct order, and finding the causative muation. I think this because none of the figures 3,4 and 9-12 show any real noticeable patterns in the data, with the exception of figure 10. The trend showed in figure 10, is not significant enough to draw a clear idea about where each fragment should be positioned, based on it's SNP skew. Therefore I conclude that the use of fragment's SNP skew, will not be a fruitful approach in rearranging fragments of the model genome, to identify the causative mutation.

For the complete R script used to generate figures 2-10 see [this commit](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/commit/5580ac551555b1541dde114e6eda1fde92ff5111) of skew_scatta.R or the [final version](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/skew_scatta.R) for figures 11/12.
