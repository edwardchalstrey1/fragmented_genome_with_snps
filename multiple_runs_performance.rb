#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'
require_relative 'lib/reform_ratio'
require 'pdist'

dataset = ARGV[0]
fasta = ReformRatio.fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")
original_order = ReformRatio.fasta_id_n_lengths(fasta)[0]

gens, fits, runs, perms = UPlot.plot_info(dataset)

UPlot.uplot(dataset, gens, fits, runs, 'big_plot1') # fitness

deviation_distances = []
perms.each do |perm|
	deviation_distances << PDist.deviation(original_order, perm)
end
UPlot.uplot(dataset, gens, deviation_distances, runs, 'big_plot_dev')

square = []
perms.each do |perm|
	square << PDist.square(original_order, perm)
end
UPlot.uplot(dataset, gens, square, runs, 'big_plot_square')

hamming = []
perms.each do |perm|
	hamming << PDist.hamming(original_order, perm)
end
UPlot.uplot(dataset, gens, hamming, runs, 'big_plot_hamming')

rdist = []
perms.each do |perm|
	rdist << PDist.rdist(original_order, perm)
end
UPlot.uplot(dataset, gens, rdist, runs, 'big_plot_r')

lcs = []
perms.each do |perm|
	lcs << PDist.lcs(original_order, perm)
end
UPlot.uplot(dataset, gens, lcs, runs, 'big_plot_lcs')

kt = []
perms.each do |perm|
	kt << PDist.kendalls_tau(original_order, perm)
end
UPlot.uplot(dataset, gens, kt, runs, 'big_plot_kendalls')