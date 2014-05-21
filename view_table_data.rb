#encoding: utf-8
require 'text-table'
require_relative 'lib/write_it'

dataset = ARGV[0]
run = ARGV[1]
generation = ARGV[2]

table_data = WriteIt.file_to_array("arabidopsis_datasets/#{dataset}/#{run}/Gen#{generation}/table_data.txt")
table_data = table_data.each_slice(4).to_a

reduced_td, population_ids = [], []
table_data.each do |permutation, fitness, type, ids|
	reduced_td << [permutation, fitness, type]
	population_ids << ids
end

puts reduced_td.to_table(:first_row_is_head => true)

best_permutation = population_ids[-1] # permutations are saved in ascending order of fitness