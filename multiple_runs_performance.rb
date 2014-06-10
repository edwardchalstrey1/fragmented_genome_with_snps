#encoding: utf-8
require_relative 'lib/score_plots/umbrella_plot'

dataset = ARGV[0]

gens, fits, runs = UPlot.plot_info(dataset)
puts 'got info'

UPlot.uplot(dataset, gens, fits, runs, 'big_plot1')
puts 'got plot'