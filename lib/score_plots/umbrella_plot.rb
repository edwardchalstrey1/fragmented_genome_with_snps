#encoding: utf-8
class UPlot
	require 'rinruby'

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

	def plot_info(dataset)
		runs = UPlot.get_dirs(dataset)
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

	# Makes plot from arrays of generations (on for each data point), metric scores, and group (the run the data is from)
	def self.uplot(generations, metric_scores, run_group)
		myr = RinRuby.new(echo = false)
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"
		myr.assign 'generations', generations
		myr.assign 'metric_scores', metric_scores
		myr.assign 'runs', run_group
		myr.eval 'uplot(generations, metric_scores, runs)'
		# myr.eval 'ggsave()'
		myr.quit
	end
end

