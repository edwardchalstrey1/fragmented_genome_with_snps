#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]
method = ARGV[1]

title = 'Many replicate runs of GATOC with varying parameter groupings'

if method == 'csv' || method == 'both'
	UPlot.data_save(dataset)
end

if method == 'plot' || method == 'both'
	shorts, x = ['fits','dev','square','ham','r_dist','lcs','kt'], 0
	['Fitness','Deviation','Square','Hamming','R','LCS','KT'].each do |metric|
		filename = "umbrella_plot_#{[shorts[x]]}"
		UPlot.uplot(dataset, filename, metric, title)
		x+=1
	end
end