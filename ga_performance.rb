#encoding: utf-8
require_relative 'score_plots/score_plots.rb'

MetricPlot::gg_plots(11, 0, 1, 'dev', 'ordinal_similarity_(deviation_distance)_gen_0-10')
MetricPlot::gg_plots(11, 0, 1, 'sq', 'square_deviation_distance_gen_0-10')
MetricPlot::gg_plots(11, 0, 1, 'ham', 'generalized_hamming_distance_gen_0-10')
MetricPlot::gg_plots(11, 0, 1, 'mod', 'modified_hamming_distance_gen_0-10')
MetricPlot::gg_plots(11, 0, 1, 'r', 'r_distance(compliment_proportion)_gen_0-10')
MetricPlot::gg_plots(11, 0, 1, 'lcs', 'longest_common_subsequence_gen_0-10')
MetricPlot::gg_plots(11, 0, 1, 'kt', 'kendalls_tau_distance_gen_0-10')