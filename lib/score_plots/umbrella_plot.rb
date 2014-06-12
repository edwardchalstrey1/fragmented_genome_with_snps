#encoding: utf-8
class UPlot
	require 'rinruby'
	require_relative 'score_plots'
	require 'pp'

	# Returns array of names of runs in the dataset
	def self.get_runs(dataset)
		runs = []
		Dir.entries("arabidopsis_datasets/#{dataset}").each do |entry|
			if entry.include?('p_run')
				runs << entry
			end
		end
		return runs
	end

	# Returns number of generations in run including 0
	def self.get_gens(dataset, run)
		gens = []
		Dir.entries("arabidopsis_datasets/#{dataset}/#{run}").each do |entry|
			if entry.include?('Gen')
				unless entry.include?('lists') || entry.include?('K')
					gens << entry
				end
			end
		end
		return gens.length
	end

	# Returns arrays of the generations, fitness scores and runs for each permutation (index in each array is from same permutation)
	# Also returns arrays of all the permutations (arrays of fragment ids), and the "parameter type", which denotes the parameters used in that run
	def self.plot_info(dataset)
		runs = UPlot.get_runs(dataset)
		gens, fitness, all_runs, all_perms, param_types = [],[],[],[],[]
		runs.each do |run|
			gen_num = UPlot.get_gens(dataset, run) # get number of generations for this run
			pops = MetricPlot.get_perms(gen_num, 0, 1, dataset, run) # all the populations of this run
			gen = 0
			pops.each do |pop|
				pop.each do |perm|
					fitness << (1 - perm[0].to_f) # compliment fitness score
					gens << gen
					all_runs << run
					all_perms << perm[1..-1]
					run_num = run.dup
					run_num.slice!('p_run')
					case run_num.to_i
					when 1..10 then param_types << 'p1'
					when 11..20 then param_types << 'p2'
					when 21..30 then param_types << 'p3'
					when 31..40 then param_types << 'p4'
					when 41..50 then param_types << 'p5'
					when 51..60 then param_types << 'p6'
					when 61..70 then param_types << 'p7'
					when 71..80 then param_types << 'p8'
					when 81..90 then param_types << 'p9'
					when 91..100 then param_types << 'p10'
					when 100..110 then param_types << 'p11'
					when 110..120 then param_types << 'p12'
					end
				end
				gen+=1
			end
		end
		return gens, fitness, all_runs, all_perms, param_types
	end

	# Makes plot from arrays of generations (on for each data point), metric scores, and the run/replicate the data is from. Runs from different parameter replicates are faceted
	def self.uplot(dataset, gens, scores, runs, param_types, filename, metric, title)
		myr = RinRuby.new(echo = false)
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"

		arrays = ['gens', 'scores', 'runs', 'param_types']
		n = 0
		[gens, scores, runs, param_types].each do |array|
			puts 'next array'
			x = 400
			myr.eval "#{arrays[n]} <- c()"
			array.each_slice(400) do |chunk|
				myr.chunk = chunk
				myr.eval "#{arrays[n]} <- c(#{arrays[n]}, chunk)"
				puts x
				x+=(chunk.length)
			end
			puts x
			n+=1
		end

		# puts 'hello'
		# myr.gens = gens
		# myr.scores = scores
		# myr.runs = runs
		# myr.param_types = param_types

		myr.dataset = dataset
		myr.filename = filename
		myr.metric = metric
		myr.title = title
		myr.eval "p <- uplot(gens, scores, runs, param_types, metric, title)"
		myr.eval "ggsave(p, file = paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset,'/', filename,'.png', sep = ''))"
		myr.quit
	end
end

