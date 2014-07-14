#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]
method = ARGV[1]

title = 'Replicate runs of a genetic algorithm, that rearranges unordered contigs from a model of a
backcrossed EMS mutagenized Arabidopsis chromosome (4). A fitness score is attributed to each permutation
of the contig order.'

unless method == nil

	if method == 'csv' || method == 'both'
		UPlot.fits_save(dataset)
	end

	if method == 'plot' || method == 'both'
		shorts = 'fits_total'
		y_axis = 'Permutation fitness score'
		filename = "umbrella_plot_#{shorts}"
		UPlot.uplot(dataset, filename, 'Fitness', y_axis, title, 'data_fits.csv')
	end

end