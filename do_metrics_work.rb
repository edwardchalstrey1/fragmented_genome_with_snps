#encoding: utf-8
require_relative 'lib/metric_work'

dataset = ARGV[0]
size = ARGV[1] # size of each population of permuations that are progressively further from correct
pop_num = ARGV[2] # number of populations
div = ARGV[3] # number of divisions of the genome at which to calculate the ratio
metric = ARGV[4]

MetricWork.adjacent_swaps_csv(dataset, size, pop_num, div, metric)
MetricWork.metric_test_plot(dataset, 'adj_test', metric, metric, 'adjacent swaps')
