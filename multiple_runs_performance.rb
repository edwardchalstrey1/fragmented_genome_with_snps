#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]
method = ARGV[1]

title = 'Replicate runs of a genetic algorithm, that rearranges unordered contigs from a model of a
backcrossed EMS mutagenized Arabidopsis chromosome (4). A fitness score is attributed to each permutation
of the contig order, based on similarity of the homozygous to heterozygous SNP ratio with an expected distribution.
Each facet contains replicates that have been run with the same parameters.'

unless method == nil

	if method == 'csv' || method == 'both'
		UPlot.data_save(dataset)
	end

	if method == 'plot' || method == 'both'
		shorts, x = ['fits','dev','square','ham','r_dist','lcs','kt'], 0
		y_axis = ['Permutation fitness score (similarity of SNP ratio across the permutation, to that of correctly ordered contigs)',
			'Similarity score (compliment proportion of deviation distance between permutations and the correct contig arrangement), normalized between 0 and 1',
			'Squared similarity score (compliment proportion of squared deviation distance between permutations and the correct contig arrangement), normalized between 0 and 1',
			'Hamming similarity (compliment proportion of hamming distance between permutations and the correct contig arrangement), normalized between 0 and 1',
			'R distance between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Longest common subsequence between permutations and the correct contig arrangement, normalized between 0 and 1',
			'Compliment of Kendall\'s tau distance between permutations and the correct contig arrangement, normalized between 0 and 1']
		['Fitness','Deviation','Square','Hamming','R','LCS','KT'].each do |metric|
			filename = "umbrella_plot_#{[shorts[x]]}"
			UPlot.uplot(dataset, filename, metric, y_axis[x], title, 'data.csv')
			x+=1
		end
	end

end