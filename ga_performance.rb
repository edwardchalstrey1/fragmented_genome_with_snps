#encoding: utf-8
require_relative 'lib/score_plots/score_plots.rb'
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

### Circos config files ###

start = s
last = s + (i*(g-1))
step = i

CircosLinks::make_links(start, last, step)
