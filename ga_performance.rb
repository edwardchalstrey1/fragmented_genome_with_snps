#encoding: utf-8
require_relative 'score_plots/score_plots.rb'
require_relative 'score_plots/example_perms.rb'
require_relative 'circos/create_links'
require 'pp'

## ARGV[0] is dataset ##
## ARGV[1] is run ##

### Plots for algorithm performance over generations ##

s = 0 # First generation in figure (start)
i = 1 # Number of generations to increment by
g = 11 # Number of generations in the plot

all_perms = MetricPlot::get_perms(g, s, i)

# MetricPlot::gg_plots(s, i, "all_metrics_gen_#{s}-#{(i*(g-1))+s}", all_perms)

MetricPlot::gg_plots(s, i, "dev", "ordinal_similarity_(deviation_distance)_gen_#{s}-#{(i*(g-1))+s}", all_perms)
MetricPlot::gg_plots(s, i, "sq", "square_deviation_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms)
MetricPlot::gg_plots(s, i, "ham", "generalized_hamming_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms)
MetricPlot::gg_plots(s, i, "mod", "modified_hamming_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms)
MetricPlot::gg_plots(s, i, "r", "r_distance(compliment_proportion)_gen_#{s}-#{(i*(g-1))+s}", all_perms)
MetricPlot::gg_plots(s, i, "lcs", "longest_common_subsequence_gen_#{s}-#{(i*(g-1))+s}", all_perms)
MetricPlot::gg_plots(s, i, "kt", "kendalls_tau_distance_gen_#{s}-#{(i*(g-1))+s}", all_perms)


### Example plots ###

# n = 10 # Number of permutations in each population (mutants of the original order, and random shuffles)

# snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf")
# fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")
# example_perms = ExamplePerms::get_perms(fasta, n, snp_data) 

# MetricPlot::gg_plots(0, 1, 'dev', '10mut_10shuf_dev', example_perms)
# MetricPlot::gg_plots(0, 1, 'sq', '10mut_10shuf_sq', example_perms)
# MetricPlot::gg_plots(0, 1, 'ham', '10mut_10shuf_ham', example_perms)
# MetricPlot::gg_plots(0, 1, 'mod', '10mut_10shuf_mod', example_perms)
# MetricPlot::gg_plots(0, 1, 'r', '10mut_10shuf_r', example_perms)
# MetricPlot::gg_plots(0, 1, 'lcs', '10mut_10shuf_lcs', example_perms)
# MetricPlot::gg_plots(0, 1, 'kt', '10mut_10shuf_kt', example_perms)


### Circos config files ###

# start = 1
# last = 3001
# step = 250

# CircosLinks::make_links(start, last, step)
