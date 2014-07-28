#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]
data_plot = ARGV[1]

title = 'Replicate runs of a genetic algorithm, that rearranges unordered contigs from a model of a
backcrossed EMS mutagenized Arabidopsis chromosome (4). A fitness score is attributed to each permutation
of the contig order.'

unless data_plot == nil

	if data_plot == 'csv' || data_plot == 'both'
		UPlot.data_save(dataset)
	end

	if data_plot == 'plot' || data_plot == 'both'
		shorts, x = ['fits','dev','square','ham','r_dist','lcs','kt'], 0
		y_axis = ['Permutation fitness score',
			'Deviation distance between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Squared deviation distance between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Hamming distance between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Compliment of R distance between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Compliment of the longest common subsequence between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Kendall\'s tau distance between permutations and the correct contig arrangement, normalized between 0 and 1']
		['Fitness','Deviation','Square','Hamming','R','LCS','KT'].each do |metric|
			filename = "umbrella_plot_#{[shorts[x]]}"
			UPlot.uplot(dataset, filename, metric, y_axis[x], title, 'data.csv', 'no_correct')
			x+=1
		end
	end

end