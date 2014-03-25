#encoding: utf-8
require_relative 'score_plots/score_plots.rb'
require_relative 'score_plots/example_perms.rb'
require 'pp'

### Plots for algorithm performance over generations ##

s = 0 # First generation in figure (start)
i = 4 # Number of generations to increment by
g = 48 # Last generation in plot

all_perms = MetricPlot::get_perms(g, s, i)

MetricPlot::gg_plots(s, i, 'dev', 'ordinal_similarity_(deviation_distance)_gen_0-48', all_perms)
MetricPlot::gg_plots(s, i, 'sq', 'square_deviation_distance_gen_0-48', all_perms)
MetricPlot::gg_plots(s, i, 'ham', 'generalized_hamming_distance_gen_0-48', all_perms)
MetricPlot::gg_plots(s, i, 'mod', 'modified_hamming_distance_gen_0-48', all_perms)
MetricPlot::gg_plots(s, i, 'r', 'r_distance(compliment_proportion)_gen_0-48', all_perms)
MetricPlot::gg_plots(s, i, 'lcs', 'longest_common_subsequence_gen_0-48', all_perms)
MetricPlot::gg_plots(s, i, 'kt', 'kendalls_tau_distance_gen_0-48', all_perms)


### Example plots ###

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf")
fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")

example_perms = ExamplePerms::get_perms(fasta, 10, snp_data) # There are 10 permutations in each population
MetricPlot::gg_plots(0, 1, 'dev', '10mut_10shuf', example_perms)


