#encoding: utf-8
require_relative 'lib/score_plots/score_plots.rb'
require_relative 'lib/score_plots/example_perms.rb'
require_relative 'circos/create_links'

location = "fragmented_genome_with_snps/arabidopsis_datasets"
dataset = ARGV[0]
run = ARGV[1]

### Plots for algorithm performance over generations ##

s = ARGV[2].to_i # First generation in figure (start)
i = ARGV[3].to_i # Number of generations to increment by
g = ARGV[4].to_i # Number of generations in the plot

all_perms = MetricPlot::get_perms(g, s, i, dataset, run)

MetricPlot::metric_plot(s, i, "dev", "ordinal_similarity_(deviation_distance)_gen_#{s}-#{(i*(g-1))+s}", all_perms, dataset, run)
MetricPlot::metric_plot(s, i, "sq", "square_deviation_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms, dataset, run)
MetricPlot::metric_plot(s, i, "ham", "generalized_hamming_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms, dataset, run)
MetricPlot::metric_plot(s, i, "r", "r_distance(compliment_proportion)_gen_#{s}-#{(i*(g-1))+s}", all_perms, dataset, run)
MetricPlot::metric_plot(s, i, "lcs", "longest_common_subsequence_gen_#{s}-#{(i*(g-1))+s}", all_perms, dataset, run)
MetricPlot::metric_plot(s, i, "kt", "kendalls_tau_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms, dataset, run)


### Example plots ###

div = ARGV[5].to_f # Number of divisions at which to calculate SNP density in permutation
n = ARGV[6].to_i # Number of permutations in each population (mutants of the original order, and random shuffles)

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{dataset}/snps.vcf")
fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")
hm = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/hm_snps.txt") # we can use the SNPs from the model genome to make example ratio
ht = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/ht_snps.txt")
fratio_breaks = SNPdist::fratio(hm, ht, div, genome_length) # frequency ratio array # 10,000 for 10K dataset
comparable_ratio = SNPdist::hyp_snps(fratio_breaks, div, genome_length) # hypothetical snp positions array

example_perms = ExamplePerms::get_perms(fasta, n, snp_data, comparable_ratio, div, genome_length) 

MetricPlot::metric_plot(0, 2, 'dev', '10_chunk_swap_shuf_dev', example_perms, dataset, run)
MetricPlot::metric_plot(0, 2, 'sq', '10_chunk_swap_shuf_sq', example_perms, dataset, run)
MetricPlot::metric_plot(0, 2, 'ham', '10_chunk_swap_shuf_ham', example_perms, dataset, run)
MetricPlot::metric_plot(0, 2, 'r', '10_chunk_swap_shuf_r', example_perms, dataset, run)
MetricPlot::metric_plot(0, 2, 'lcs', '10_chunk_swap_shuf_lcs', example_perms, dataset, run)
MetricPlot::metric_plot(0, 2, 'kt', '10_chunk_swap_shuf_kt', example_perms, dataset, run)


### Circos config files ###

# start = 1
# last = 3001
# step = 250

# CircosLinks::make_links(start, last, step)
