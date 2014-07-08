#encoding: utf-8
require_relative 'lib/score_plots/score_plots.rb'
require_relative 'circos/create_links'

### Make Circos config files ###

location = "fragmented_genome_with_snps/arabidopsis_datasets"
dataset = ARGV[0]
run = ARGV[1]
s = ARGV[2].to_i # First generation in figure (start)
i = ARGV[3].to_i # Number of generations to increment by
g = ARGV[4].to_i # Number of generations in the plot

last = s + (i*(g-1))

CircosLinks::make_links(s, last, i)