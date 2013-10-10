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
def odds (array) #create array from odd indexes of input array
	return array.values_at(*array.each_index.select {|i| i.odd?})
end
def evens (array)
	return array.values_at(*array.each_index.select {|i| i.even?})
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

# Even/odd and odd/even
def even_odd (frags_by_density, odd_even_string)
	if odd_even_string == 'even' #even values then odd values reversed or opposite
		rearranged = evens(frags_by_density) #the even numbers first half
		rearranged << (odds(frags_by_density)).reverse #the odd numbers reversed second half 
	else
		rearranged = odds(frags_by_density) #the odd numbers first half
		rearranged << (evens(frags_by_density)).reverse #the even numbers reversed second half 
	end
	return rearranged.flatten   
end

# Left/right method
def left_right (id_pos_hash, fasta_lengths)
	left = [] # frags with snps's mostly to the right, so should be on left of distribution
	right = [] # frags with snps's mostly to the left, so should be on right of distribution
	d = 0
	id_pos_hash.keys.each do |x| # for each of the frag ids
		temp = []
		(id_pos_hash[x].split(",").map { |s| s.to_i }).each do |i| # adding each of the positions associated with that frag/key to temporary array
			temp << i
		end
		sum = temp.inject(:+) # summing the positions
		if sum != nil # ignoring the fragments with 0 snps, which have no postions, so sum = nil
			lxnp = (fasta_lengths[d]*temp.length)/2 # then working out half the length of the fragment, multiplied by the number of positions
			if sum < lxnp # if the sum of the positions is < half the length of the fragment, multiplied by the number of positions
				# this is equivalent to working out the average position and checking it against half the length
				right << x
			else
				left << x
			end
		else #frags with 0 snps
			if evens(id_pos_hash.keys).include?(x) 
				left << x
			else
				right << x #add odd indexed frags with 0 snps to right, and evens to left
			end
		end 
		d+=1 
	end
	return left, right
end




frags_original_order = extract_json('frag_ids_original_order.json')
frags_by_density = extract_json('frags_by_density.json')

id_pos_hash = extract_json('pos_hash.json') #remember the fragments are in a random order here - the order they were in in the fasta
fasta_lengths = extract_json('fasta_lengths.json') #also in the order from the fasta file

#puts 'Density order Score: ' + (score(frags_original_order, frags_by_density)).to_s
#puts 'Random Score: ' + (random_score(frags_original_order, frags_by_density)).to_s
puts 'Even Odd Method Score: ' + (score(frags_original_order, (even_odd(frags_by_density, 'even')))).to_s
puts 'Odd Even Method Score: ' + (score(frags_original_order, (even_odd(frags_by_density, 'odd')))).to_s
puts 'Left Right Method Score: ' + (score(frags_original_order, (left_right(id_pos_hash, fasta_lengths).flatten))).to_s