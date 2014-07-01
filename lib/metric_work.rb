#encoding: utf-8

class MetricWork

	require_relative 'reform_ratio'
	require_relative 'fitness_score'
	require_relative 'GATOC'
	require 'pmeth'
	require 'pdist'

=begin Testing how well the fitness method of the genetic algorithm,
		identifies permuations that are approaching the correct order.

		1. Get an array of correctly ordered FASTA contigs
		2. Get the ratio of the correctly ordered contigs, so we can workin out fitness of permutation
		3. Create permutations that are progressively further from the correct
		4. Save the permutation fitness scores directly into a csv, to make a ggplot (but can also save permutation txt files - possibly not neccesary)
		5. Create a plot to represent this
=end

	def self.score(fasta, perm, snp_data, comparable_ratio, div, genome_length, metric)
		case metric
		when 'Fitness'
			score = GATOC.fitness(perm, snp_data, comparable_ratio, div, genome_length)[0]
		when 'DeviationDistance'
			score = PDist.deviation(fasta, perm)
		when 'SquareDeviationDistance'
			score = PDist.square(fasta, perm)
		when 'HammingDistance'
			score = PDist.hamming(fasta, perm)
		when 'RDistance'
			score = PDist.rdist(fasta, perm)
		when 'LongestCommonSubsequence'
			score = PDist.lcs(fasta, perm)
		when 'KendallsTau'
			score = PDist.kendalls_tau(fasta, perm)
		end
		return score
	end

	def self.adjacent_swaps_csv(dataset, size, pop_num, div, metric, filename)
		# 1
		fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta") # correct permutation
		###

		# 2
		snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{dataset}/snps.vcf")
		genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")
		ht, hm = ReformRatio.perm_pos(fasta, snp_data)
		comparable_ratio = FitnessScore::ratio(hm, ht, div, genome_length)
		###

		# 3/4: Population of adjacent_swap mutants (of the correct contig order)
		start_pop = []
		size.times do
		start_pop << fasta
		end

		WriteIt.add_to("arabidopsis_datasets/#{dataset}/#{filename}.csv", "population,#{metric}")
		x = 1
		pop_num.times do
			adj_pop = []
			start_pop.each do |perm|
				new_perm = PMeth.adjacent_swap(perm)
				adj_pop << new_perm # need this population to be the next starting population
				score = MetricWork.score(fasta, perm, snp_data, comparable_ratio, div, genome_length, metric)
				WriteIt.add_to("arabidopsis_datasets/#{dataset}/#{filename}.csv", "#{x},#{score}")
			end
			start_pop = adj_pop
			x+=1
		end
	end

	# 5:
	def self.metric_test_plot(dataset, filename, metric, x_axis, y_axis, title, input_file)
		myr = RinRuby.new(echo = false)
		myr.dataset = dataset
		myr.filename = filename
		myr.metric = metric
		myr.title = title
		myr.x_axis = x_axis
		myr.y_axis = y_axis
		myr.input_file = input_file
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"
		myr.eval "df <- read.csv(paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset, '/', input_file, '.csv', sep=''))"
		myr.eval "p <- metric_test_plot(df, title, x_axis, y_axis, metric)"
		myr.eval "ggsave(p, file = paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset,'/', filename,'.png', sep = ''))"
		myr.quit
		puts 'made a plot'
	end
end
