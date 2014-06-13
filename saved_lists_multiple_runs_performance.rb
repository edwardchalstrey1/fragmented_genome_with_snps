#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]
title = 'Many replicate runs of GATOC with varying parameter groupings'

gens = WriteIt.file_to_ints_array("arabidopsis_datasets/#{dataset}/gens.txt")
runs = WriteIt.file_to_array("arabidopsis_datasets/#{dataset}/runs.txt")
param_types = WriteIt.file_to_array("arabidopsis_datasets/#{dataset}/param_types.txt")
fits = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/fits.txt")
deviation_distances = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/dev.txt")
square = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/sq.txt")
hamming = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/ham.txt")
rdist = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/rdist.txt")
lcs = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/lcs.txt")
kt = WriteIt.file_to_floats_array("arabidopsis_datasets/#{dataset}/kt.txt")

UPlot.uplot(dataset, gens, fits, runs, param_types, 'big_plot_fit', 'Fitness (compliment)', title)
UPlot.uplot(dataset, gens, deviation_distances, runs, param_types, 'big_plot_dev', 'Deviation distance (normalized)', title)
UPlot.uplot(dataset, gens, square, runs, param_types, 'big_plot_square', 'Square deviation distance (normalized)', title)
UPlot.uplot(dataset, gens, hamming, runs, param_types, 'big_plot_hamming', 'Hamming distance (normalized)', title)
UPlot.uplot(dataset, gens, rdist, runs, param_types, 'big_plot_r', 'R distance (normalized compliment)', title)
UPlot.uplot(dataset, gens, lcs, runs, param_types, 'big_plot_lcs', 'Longest common subsequence (normalized)', title)
UPlot.uplot(dataset, gens, kt, runs, param_types, 'big_plot_kendalls', 'Kendall\'s tau (normalized)', title)