#encoding: utf-8
require_relative 'lib/score_plots/score_plots'
require_relative 'lib/score_plots/umbrella_plot'

dataset = '10K_dataset3'

runs = UPlot.get_dirs(dataset)

gens, fitness, runs = [],[],[]

runs.each do |run|
	gen_num = UPlot.get_gens(dataset, run) # get number of generations for this run
	pops = MetricPlot.get_perms(gen_num, 0, 1, dataset, run) # all the populations of this run
	gen = 0
	pops.each do |pop|
		pop.each do |perm|
			fitness << perm[0]
			gens << gen
			runs << 
		end
		gen+=1
	end