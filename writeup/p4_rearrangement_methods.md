Rearrangement methods part 4
========================================================

### Testing the SNP Skew

Since the Kernel density estimate was not normalized for each fragment, I have come up with a new way of testing SNP skew as a method for rearranging fragments with SNPs into the correct order. See the project [README](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/README.md) and parts [1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/rearrangement_methods.md), [2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/p2_rearrangement_methods.md) and [3](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/writeup/p3_rearrangement_methods.md) of rearrangement methods for details.

The new method I have used for determining SNP gradients for each of the fragments, relies on a different way of representing SNP density across the fragment, which I will refer to as the distance method. The distance method...

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_abs_2_10000.png?raw=true)
[Figure 1](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_abs_2_10000.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_grad_2_10000.png?raw=true)
[Figure 2](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_grad_2_10000.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/example_gradient_f681_mins2.png?raw=true)
[Figure 3](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/example_gradient_f681_mins2.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/example_gradient_f258_mins2.png?raw=true)
[Figure 4](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/example_gradient_f258_mins2.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_abs_30_10000.png?raw=true)
[Figure 5](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_abs_30_10000.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_grad_30_10000.png?raw=true)
[Figure 6](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/skew_scatter_grad_30_10000.png)

![Image](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/example_gradient_f729_mins30.png?raw=true)
[Figure 7](https://github.com/edwardchalstrey1/fragmented_genome_with_snps/blob/master/arabidopsis_datasets/dataset5/figures/example_gradient_f729_mins30.png)
