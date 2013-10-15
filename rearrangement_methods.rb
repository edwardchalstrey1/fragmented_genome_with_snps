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
def two_a_to_h (keys_array, values_array)
	return Hash[*keys_array.zip(values_array).flatten]
end

def density_order (frags_by_density, fasta_lengths, id_pos_hash) #get fasta lengths array and id/positions hash in SNP density order
	d_lengths = []
	d_pos = []
	length_hash = two_a_to_h(id_pos_hash.keys, fasta_lengths) #make hash of ids in fasta file order with lengths
	frags_by_density.each do |id|
		d_pos << id_pos_hash[id]
		d_lengths << length_hash[id] #the lengths of each fragment in the same order as frags_by_density (SNP density order)
	end
	d_id_pos_hash = two_a_to_h(frags_by_density, d_pos) #an id/positions hash where the ids are in SNP density order
	return d_id_pos_hash, d_lengths
end

# Rearrangement Methods
# ---------------------

# Control
def random_score (frags_original_order, frags_by_density)
	scores = []
	100.times do |i|
		scores << score(frags_original_order, frags_by_density.shuffle)
	end
	av_score = scores.inject(:+).to_f / 100
	return scores
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
def lr_d (d_id_pos_hash, d_lengths)
	lnr = left_right(d_id_pos_hash, d_lengths) #get the left/right arrays for density ordered frags so both arrays are still ordered by density
	left = lnr[0]
	right = lnr[1]
	left << right.reverse # reverse the right array into descending density order
	return left.flatten
end
def max_score (original_order)
end
def re_order_densities (id_density_hash, rearranged_ids) # get list of the densities for a rearrangement of ids (frags)
	densities = []
	rearranged_ids.each do |id|
		densities << id_density_hash[id]
	end
	return densities
end


frags_original_order = extract_json('frag_ids_original_order.json')
frags_reverse_order = frags_original_order.reverse
frags_by_density = extract_json('frags_by_density.json')

id_pos_hash = extract_json('pos_hash.json') #remember the fragments are in a random order here - the order they were in in the fasta
fasta_lengths = extract_json('fasta_lengths.json') #also in the order from the fasta file

density_order_data = density_order(frags_by_density, fasta_lengths, id_pos_hash)
d_id_pos_hash = density_order_data[0]
d_lengths = density_order_data[1]

id_density_hash = extract_json('id_density_hash.json')

#puts 'Density order Score: ' + (score(frags_original_order, frags_by_density)).to_s
#puts 'Random Score: ' + (random_score(frags_original_order, frags_by_density)).to_s
#puts 'Even Odd Method Score: ' + (score(frags_original_order, (even_odd(frags_by_density, 'even')))).to_s
#puts 'Odd Even Method Score: ' + (score(frags_original_order, (even_odd(frags_by_density, 'odd')))).to_s
#puts 'Left Right Method Score: ' + (score(frags_original_order, (left_right(id_pos_hash, fasta_lengths).flatten))).to_s
#puts 'Left Right Density Method Score: ' + (score(frags_original_order, (lr_d(d_id_pos_hash, d_lengths)))).to_s
#puts score(frags_reverse_order, frags_original_order)


File.open('random_scores.txt', "w+") do |f|
	(random_score(frags_original_order, frags_by_density)).each { |i| f.puts(i) }
end
