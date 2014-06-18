#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]

title = 'Many replicate runs of GATOC with varying parameter groupings'

UPlot.data_save(dataset)
# ['fitness', 'dev', 'square', 'ham', 'r_dist', 'lcs', 'kt'].each do |metric|
# 	filename = "umbrella_plot_#{metric}"
# 	UPlot.uplot(dataset, filename, metric, title)
# end