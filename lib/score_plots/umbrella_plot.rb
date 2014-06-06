#encoding: utf-8
class UPlot
	require 'rinruby'

	def self.get_dirs(dataset)
		run_dirs = [] # these are the directories we will be taking the information from
		Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}")) do
			dirs = Dir.glob('*').select {|f| File.directory? f}
			dirs.each do |dir|
				if dir.include?('p_run')
					run_dirs << "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}/#{dir}"
				end
			end
		end
		return run_dirs
	end

	# def self.

	def self.uplot(generations, metric_scores, runs)
		myr = RinRuby.new(echo = false)
		myr.eval "source('~/fragmented_genome_with_snps/lib/score_plots/umbrella_plot.R')"
		myr.assign 'generations', generations
		myr.assign 'metric_scores', metric_scores
		myr.assign 'runs', runs
		myr.eval 'uplot(generations, metric_scores, runs)'
		# myr.eval 'ggsave()'
		myr.quit
	end
end

