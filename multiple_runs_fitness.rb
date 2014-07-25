#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'
require_relative 'lib/GATOC'
require_relative 'lib/reform_ratio'

dataset = ARGV[0]
data_plot = ARGV[1]
fitness_method = ARGV[2]

title = 'Replicate runs of a genetic algorithm, that rearranges unordered contigs from a model of a
backcrossed EMS mutagenized Arabidopsis chromosome (4). A fitness score is attributed to each permutation
of the contig order.'

if fitness_method != nil
	fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
	fasta = ReformRatio.fasta_array(fasta_file)
	snp_data = ReformRatio.get_snp_data("arabidopsis_datasets/#{dataset}/snps.vcf")
	genome_length = ReformRatio.genome_length(fasta_file)
	correct_fitness = GATOC.fitness(fasta, snp_data, genome_length, fitness_method)[0]
else
	correct_fitness = 'no_correct'
end

unless data_plot == nil

	if data_plot == 'csv' || data_plot == 'both'
		UPlot.fits_save(dataset)
	end

	if data_plot == 'plot' || data_plot == 'both'
		shorts = 'fits_total'
		y_axis = 'Permutation fitness score'
		if fitness_method == nil
			filename = "umbrella_plot_#{shorts}"
		else
			filename = "umbrella_plot_#{shorts}_with_correct"
		end
		UPlot.uplot(dataset, filename, 'Fitness', y_axis, title, 'data_fits.csv', correct_fitness)
	end

end