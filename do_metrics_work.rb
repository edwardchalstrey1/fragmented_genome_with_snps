#encoding: utf-8
require_relative 'lib/metric_work'

dataset = ARGV[0]
div = ARGV[1].to_f # number of divisions of the genome at which to calculate the ratio
size = ARGV[2].to_i # size of each population of permuations that are progressively further from correct
pop_num = ARGV[3].to_i # number of populations
metric = ARGV[4]

# MetricWork.adjacent_swaps_csv(dataset, size, pop_num, div, metric)
MetricWork.metric_test_plot(dataset, 'adj_test', metric, metric, 'adjacent swaps')
