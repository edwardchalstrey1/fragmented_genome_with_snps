#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'
require_relative 'lib/reform_ratio'
require 'pdist'

dataset = ARGV[0]
fasta = ReformRatio.fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")
original_order = ReformRatio.fasta_id_n_lengths(fasta)[0]

title = 'Many replicate runs of GATOC with varying parameter groupings'

gens, fits, runs, perms, param_types = UPlot.plot_info(dataset)
UPlot.uplot(dataset, gens, fits, runs, param_types, 'big_plot_fit', 'Fitness (compliment)', title)

deviation_distances, x = [], 1
perms.each do |perm|
	deviation_distances << PDist.deviation(original_order, perm)
	puts "Permutations: #{x}"; x+=1
end
UPlot.uplot(dataset, gens, deviation_distances, runs, param_types, 'big_plot_dev', 'Deviation distance (normalized)', title)

square, x = [], 1
perms.each do |perm|
	square << PDist.square(original_order, perm)
	puts "Permutations: #{x}"; x+=1
end
UPlot.uplot(dataset, gens, square, runs, param_types, 'big_plot_square', 'Square deviation distance (normalized)', title)

hamming, x = [], 1
perms.each do |perm|
	hamming << PDist.hamming(original_order, perm)
	puts "Permutations: #{x}"; x+=1
end
UPlot.uplot(dataset, gens, hamming, runs, param_types, 'big_plot_hamming', 'Hamming distance (normalized)', title)

rdist, x = [], 1
perms.each do |perm|
	rdist << PDist.rdist(original_order, perm)
	puts "Permutations: #{x}"; x+=1
end
UPlot.uplot(dataset, gens, rdist, runs, param_types, 'big_plot_r', 'R distance (normalized compliment)', title)

lcs, x = [], 1
perms.each do |perm|
	lcs << PDist.lcs(original_order, perm)
	puts "Permutations: #{x}"; x+=1
end
UPlot.uplot(dataset, gens, lcs, runs, param_types, 'big_plot_lcs', 'Longest common subsequence (normalized)', title)

kt, x = [], 1
perms.each do |perm|
	kt << PDist.kendalls_tau(original_order, perm)
	puts "Permutations: #{x}"; x+=1
end
UPlot.uplot(dataset, gens, kt, runs, param_types, 'big_plot_kendalls', 'Kendall\'s tau (normalized)', title)

WriteIt.write_txt("arabidopsis_datasets/#{dataset}/gens", gens)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/runs", runs)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/param_types", param_types)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/fits", fits)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/dev", deviation_distances)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/sq", square)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/ham", hamming)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/rdist", rdist)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/lcs", lcs)
WriteIt.write_txt("arabidopsis_datasets/#{dataset}/kt", kt)