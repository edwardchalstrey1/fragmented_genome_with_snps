#encoding: utf-8
class MetricPlot

	require_relative '../write_it'
	require 'pdist'
	require_relative '../reform_ratio'
	require 'rinruby'

	# Inpput: Dataset
	# Output: The correctly ordered permutation of frag ids
	def self.original_order(dataset)
		original_fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")
		return ReformRatio::fasta_id_n_lengths(original_fasta)[0]
	end

	# Input 0: Number of generations to plot
	# Input 1: The generation to start with
	# Input 2: The number of generations to increment by
	# Input 3: Dataset
	# Input 4: Run
	# Output: Array of populations, each population is from the next generation and has a number of permutations (arrays of fasta frag ids)
	def self.get_perms(gen, start, inc, dataset, run) # number of generations to choose from
		all_perms = []
		n = start
		gen.times do
			pop = []
			Dir.entries("arabidopsis_datasets/#{dataset}/#{run}/Gen#{n}").each do |ptxt|
				unless ptxt.include?('best') || ptxt.include?('table') # excluding table_data.txt and best_permutation.txt
					if ptxt.include?('.txt')
						perm = []
						IO.foreach("arabidopsis_datasets/#{dataset}/#{run}/Gen#{n}/#{ptxt}") { |line| perm << line.gsub(/\n/,'') }
						pop << perm
					end
				end
			end
			all_perms << pop
			n+=inc
		end
		return all_perms
	end

	# Input 0: The generation to start with
	# Input 1: The number of generations to increment by
	# Input 2: String with the id of the metric to plot (see case below)
	# Input 3: Array of populations, each population is from the next generation and has a number of permutations (arrays of fasta frag ids)
	# Input 4: Dataset
	# Input 5: Run
	# Output 0: Array of generation numbers for each population
	# Output 1: Array of average metric scores for each population (>1 per population due to >1 metrics)
	# Output 2: Array of standard errors for each metric score
	# Output 3: Array of strings indicating what metric each score is from
	# Output 4: Array of best scores for each population (the fitness score is taken from the permutation with the highest metric score)
	def self.plot_info(start, inc, met, all_perms, dataset, run)
		orig = original_order(dataset)
		myr = RinRuby.new(echo = false)
		myr.eval 'source(paste(getwd(),"/lib/score_plots/score_plots.R",sep=""))'
		x, y, se, best_sc = [], [], [], []
		n = start
		all_perms.each do |pop|
			fitness, metric = [], []
			pop.each do |perm|
				fitness << (1.0 - perm[0].to_f) # the compliment proportion of the fitness value
				red_perm = perm[1..-1] # red_perm is just the permutation, reduced in length by one to not have the fitness float
				case met
				when 'dev'
					metric << PDist.deviation(orig, red_perm)
				when 'sq'
					metric << PDist.square(orig, red_perm)
				when 'ham'
					metric << PDist.hamming(orig, red_perm)
				when 'r'
					metric << PDist.rdist(orig, red_perm)
				when 'lcs'
					metric << PDist.lcs(orig, red_perm)
				when 'kt'
					metric << PDist.kendalls_tau(orig, red_perm)
				end
			end
			pop_y = [fitness, metric]
			pop_y.each do |scores|
				y << scores.inject(:+) / scores.length.to_f
				x << n
				myr.assign 'scores', scores
				sem = myr.pull 'st_err(scores)'
				se << sem.to_f
			end
			best_sc << fitness.sort[0].to_f # The best fitness score, lowest is best
			best_sc << metric[fitness.index(fitness.sort[0])].to_f # The metric score of the same permutation that has best fitness
			n+=inc
		end
		group = ['comp_fit', met] * all_perms.length
		myr.quit
		return x, y, se, group, best_sc
	end

	# Input: Inputs for plot_info, plus name of the file string at input 3, moving subsequent areguments one place along
	# Output: ggplot of the fitness over generations of the genetic algorithm, with a permutation distance metric for comparison
	def self.metric_plot(start, inc, met, filename, all_perms, dataset, run)
		x, y, se, group, best_sc = plot_info(start, inc, met, all_perms, dataset, run)
		myr = RinRuby.new(echo = false)
		myr.eval 'source(paste(getwd(),"/lib/score_plots/score_plots.R",sep=""))'
		myr.assign 'x', x
		myr.assign 'y', y
		myr.assign 'se', se
		myr.assign 'group', group
		myr.assign 'best_sc', best_sc
		myr.assign 'dataset_run', "#{dataset}/#{run}"
		myr.assign 'filename', filename
		myr.eval 'plot_it(x, y, se, group, best_sc, dataset_run, filename)'
		myr.quit
	end
end