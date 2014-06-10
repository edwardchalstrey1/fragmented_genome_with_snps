#encoding: utf-8
class UPlot
	require 'rinruby'
	require_relative 'score_plots'
	require 'pdist'
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
	def self.plot_info(dataset)
		runs = UPlot.get_runs(dataset)
		gens, fitness, all_runs = [],[],[]
		runs.each do |run|
			gen_num = UPlot.get_gens(dataset, run) # get number of generations for this run
			pops = MetricPlot.get_perms(gen_num, 0, 1, dataset, run) # all the populations of this run
			gen = 0
			pops.each do |pop|
				pop.each do |perm|
					fitness << perm[0]
					gens << gen
					all_runs << run
				end
				gen+=1
			end
		end
		return gens, fitness, all_runs
	end

	def self.get_metric

	end

	# Makes plot from arrays of generations (on for each data point), metric scores, and group (the run the data is from)
	def self.uplot(dataset, gens, scores, runs, filename)
		myr = RinRuby.new(echo = false)
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"
		puts 'a'
		arrays = ['gens', 'scores', 'runs']
		n = 0
		[gens, scores, runs].each do |array|
			x = 1
			myr.eval "#{arrays[n]} <- c()"
			array.each do |entry|
				if n == 1
					myr.assign 'entry', entry.to_f
				else
					myr.assign 'entry', entry
				end
				myr.eval "#{arrays[n]} <- c(#{arrays[n]}, entry)"
				puts x
				x+=1
			end
			n+=1
		end
		myr.assign 'dataset', dataset
		myr.assign 'filename', filename
		myr.eval "p <- uplot(gens, scores, runs)"
		puts 'got p'
		myr.eval "ggsave(p, file = paste('~/fragmented_genome_with_snps/arabidopsis_datasets/', dataset,'/', filename,'.png', sep = ''))"
		myr.quit
	end
end

