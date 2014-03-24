#encoding: utf-8
class MetricPlot

	require '~/fragmented_genome_with_snps/lib/write_it'
	require '~/fragmented_genome_with_snps/lib/rearrangement_score'
	require 'rinruby'

	# Output: The correctly ordered permutation of frag ids
	def self.original_order
		original_order = []
		IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/correct_permutation.txt") { |line| original_order << line }
		original_order[1..-1]
	end

	# Input 0: Number of generations to plot
	# Input 1: The generation to start with
	# Input 2: The number of generations to increment by
	# Output: Array of populations, each population is from the next generation and has a number of permutations (arrays of fasta frag ids)
	def self.get_perms(gen, start, inc) # number of generations to choose from
		all_perms = []
		n = start
		gen.times do
			pop = []
			Dir.entries("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{n}").each do |ptxt|
				if ptxt.include? '.txt'
					perm = []
					IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{n}/#{ptxt}") { |line| perm << line }
					pop << perm
				end
			end
			all_perms << pop
			n+=inc
		end
		return all_perms
	end

	# Input 0: Number of generations to plot
	# Input 1: The generation to start with
	# Input 2: The number of generations to increment by
	# Input 3: String with the id of the metric to plot (see case below)
	# Input 4: String - name of the plot
	# Input 5: Array of populations, each population is from the next generation and has a number of permutations (arrays of fasta frag ids)
	# Input 6: Type: example or performance figures
	# Output: ggplot of the fitness over generations of the genetic algorithm, with a permutation distance metric for comparison
	def self.gg_plots(gen, start, inc, met, filename, all_perms, type)
		orig = original_order
		# all_perms = get_perms(gen, start, inc)
		myr = RinRuby.new(echo = false)
		myr.eval 'source("~/fragmented_genome_with_snps/score_plots/score_plots.R")'
		x, y, se = [], [], []
		n = start
		all_perms.each do |pop|
			fitness, metric = [], []
			pop.each do |perm|
				if type != 'example'
					fitness << ### TODO this needs to be fitness of example
				else
					fitness << (1 - (perm[0].gsub(/\n/, "")).to_f) # the compliment proportion of the fitness value
					perm = perm[1..-1]
				end
				case met
				when 'dev'
					metric << RearrangementScore::dev_dist(orig, perm) ### TODO make sure what is going in here makes sense for example
				when 'sq'
					metric << RearrangementScore::sq_dev_dist(orig, perm)
				when 'ham'
					metric << RearrangementScore::gen_ham_dist(orig, perm)
				when 'mod'
					metric << RearrangementScore::mod_ham_dist(orig, perm)
				when 'r'
					metric << RearrangementScore::r_dist(orig, perm)
				when 'lcs'
					metric << RearrangementScore::lcs(orig, perm)
				when 'kt'
					metric << RearrangementScore::kendalls_tau(orig, perm)
				end
			end
			pop_y = [fitness, metric]
			pop_y.each do |scores|
				y << scores.inject(:+) / scores.length.to_f
				x << n
				myr.assign 'scores', scores
				sem = myr.pull 'st_err(scores)'
				se << sem
			end
			n+=inc
		end
		group = ['comp_fit', met] * all_perms.length
		myr.assign 'x', x
		myr.assign 'y', y
		myr.assign 'se', se
		myr.assign 'group', group
		myr.assign 'dataset_run', "#{ARGV[0]}/#{ARGV[1]}"
		myr.assign 'filename', filename
		myr.eval 'plot_it(x, y, se, group, dataset_run, filename)'
		myr.quit
	end
end