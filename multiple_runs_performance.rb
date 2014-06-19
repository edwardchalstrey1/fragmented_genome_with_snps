#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]

title = 'Many replicate runs of GATOC with varying parameter groupings'

UPlot.data_save(dataset)

# shorts, x = ['fits','dev','square','ham','r_dist','lcs','kt'], 0
# ['Fitness','Deviation','Square','Hamming','R','LCS','KT'].each do |metric|
# 	filename = "umbrella_plot_#{[shorts[x]]}"
# 	UPlot.uplot(dataset, filename, metric, title)
# 	x+=1
# end