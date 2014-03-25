#encoding: utf-8
require_relative 'score_plots/score_plots.rb'
require_relative 'score_plots/example_perms.rb'
require 'pp'

### Example plots ###

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf")
fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")

example_perms = ExamplePerms::get_perms(fasta, 2, snp_data)
MetricPlot::gg_plots(0, 1, 'dev', 'mut_shuf_dev', example_perms)


### Plots for algorithm performance over generations ##

# all_perms = MetricPlot::get_perms(11, 0, 1)

# MetricPlot::gg_plots(0, 1, 'dev', 'ordinal_similarity_(deviation_distance)_gen_0-10', all_perms)
# MetricPlot::gg_plots(0, 1, 'sq', 'square_deviation_distance_gen_0-10', all_perms)
# MetricPlot::gg_plots(0, 1, 'ham', 'generalized_hamming_distance_gen_0-10', all_perms)
# MetricPlot::gg_plots(0, 1, 'mod', 'modified_hamming_distance_gen_0-10', all_perms)
# MetricPlot::gg_plots(0, 1, 'r', 'r_distance(compliment_proportion)_gen_0-10', all_perms)
# MetricPlot::gg_plots(0, 1, 'lcs', 'longest_common_subsequence_gen_0-10', all_perms)
# MetricPlot::gg_plots(0, 1, 'kt', 'kendalls_tau_distance_gen_0-10', all_perms)




