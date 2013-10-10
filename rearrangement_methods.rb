require 'rubygems'
require 'bio-samtools'
require 'bio'
require "json"

def extract_json (json)
	output_array = JSON.parse(File.open(json).read)
	return output_array
end
def score (frags_original_order, rearranged)
	position_each_frag_id_in_d = frags_original_order.map{|x| rearranged.index(x)} #works out the index of fasta_id values in frags_by_density
	index_values = Array(0..(frags_original_order.length - 1)) # index values that fasta_ids originally at
	both = []
	both << position_each_frag_id_in_d
	both << index_values
	difference = both.transpose.map {|x| x.reduce(:-)} # taking away old position from new position, to find the distance that the frag has moved when re-ordered
	difference_abs = []
	difference.each do |i|
		difference_abs << i.abs
	end
	score = difference_abs.inject(:+) #high score = bad, score of 0 means the fragments in the right order
	return score
end

# Rearrangement Methods
# ---------------------

# Control
def random_score (frags_original_order, frags_by_density)
	scores = []
	100.times do |i|
		scores << score(frags_original_order, frags_by_density.shuffle)
	end
	return scores.inject(:+).to_f / 100
end

frags_original_order = extract_json('frag_ids_original_order.json')
frags_by_density = extract_json('frags_by_density.json')
id_pos_hash = extract_json('pos_hash.json') #remember the fragments are in a random order here

puts 'Density order Score: ' + (score(frags_original_order, frags_by_density)).to_s
puts 'Random Score: ' + (random_score(frags_original_order, frags_by_density)).to_s
