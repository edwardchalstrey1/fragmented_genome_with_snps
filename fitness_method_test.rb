#encoding: utf-8
require_relative 'lib/reform_ratio'
require 'pmeth'

=begin Testing how well the fitness method of the genetic algorithm,
		identifies permuations that are approaching the correct order.

		1. Get an array of correctly ordered FASTA contigs
		2. Create permutations that are progressively further from the correct
		3. Save the permutation fitness scores directly into a csv, to make a ggplot (but can also save permutation txt files - possibly not neccesary)
		4. Create a plot to represent this
=end

dataset = ARGV[0]
size = ARGV[1] # size of each population of permuations that are progressively further from correct
pop_num = ARGV[2] # number of populations 

# 1
fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta") # correct permutation
###

# 2/3: Population of adjacent_swap mutants (of the correct contig order)
start_pop = []
size.times do
	start_pop << fasta
end

pop_num.times do
	adj_pop = []
	start_pop.each do |perm|
		new_perm = PMeth.adjacent_swap(perm)
		adj_pop << new_perm
		# TODO get the permutation's fitness
		# TODO add the fitness to a csv
	end
	start_pop = adj_pop
end

# 4:
