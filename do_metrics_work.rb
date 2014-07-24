#encoding: utf-8
require_relative 'lib/metric_work'
require_relative 'lib/reform_ratio'
require_relative 'lib/fitness_score'

# Required:
dataset = ARGV[0]
size = ARGV[1].to_i # size of each population of permuations that are progressively further from correct
pop_num = ARGV[2].to_i # number of populations
data_plot = ARGV[3]
metric = ARGV[4]
swap_num = ARGV[5].to_i # number on adjacent swaps performed on permutations between each population

# Optional
div = ARGV[6].to_i # For count_ratio, max_hyp and hyp_distances fitness methods
if metric == 'CountRatio' # For count_ratio fitness method
	vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
	fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
	fasta = ReformRatio::fasta_array(fasta_file) # array of fasta format fragments
	snp_data = ReformRatio.get_snp_data(vcf_file)
	ht, hm = ReformRatio.perm_pos(fasta, snp_data)
	genome_length = ReformRatio::genome_length(fasta_file)
	expected_ratios = FitnessScore.ratio(hm, ht, div, genome_length) # Array of expected ratios (floats) of homozygous to heterozygous SNPs for each division of the genome
else
	expected_ratios = nil
end
params = {:div => div, :expected_ratios => expected_ratios}

# Plotting
x_axis = "Number of adjacent swaps carried out on each permutation"
y_axis = metric.gsub(/(?<=[A-Za-z])(?=[A-Z])/, ' ')
filename = "adjacent_swaps_#{metric}_#{pop_num}pop_#{size}size_swap#{swap_num}"
title = "The change in #{y_axis} of permuations, increasingly distant
		in the search space from the optimal arrangement.
		Permutations are of genomic contig order. Averages are taken
		for populations of permuations with the same number of adjacent
		swap mutations"

if data_plot == 'csv' || data_plot == 'both'
	MetricWork.adjacent_swaps_csv(dataset, size, pop_num, metric, filename, swap_num, params)
end
if data_plot == 'plot' || data_plot == 'both'
	MetricWork.metric_test_plot(dataset, filename, metric, x_axis, y_axis, title, filename)
end
