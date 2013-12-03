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
def re_order_densities (id_density_hash, rearranged_ids) # get list of the densities for a rearrangement of ids (frags)
	densities = []
	rearranged_ids.each do |id|
		densities << id_density_hash[id]
	end
	return densities
end
def scatta_txts (array, name_string) # write each list of ids and densities into id_n_density_txts directory, for use in scatter plots in R
	File.open(('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/id_n_density_txts/'+name_string+'.txt'), "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end
def write_txt (filename, array)
	File.open(filename, "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end
# Rearrangement Methods
# ---------------------

# Control
def random_score (frags_original_order, frags_by_density)
	scores = []
	100.times do |i|
		scores << score(frags_original_order, frags_by_density.shuffle)
	end
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
def lr_d (frags_by_density, fasta_lengths, id_pos_hash)
	density_order_data = density_order(frags_by_density, fasta_lengths, id_pos_hash)
	d_id_pos_hash = density_order_data[0]
	d_lengths = density_order_data[1]
	lnr = left_right(d_id_pos_hash, d_lengths) #get the left/right arrays for density ordered frags so both arrays are still ordered by density
	return lnr[0], lnr[1].reverse # reverse the right array into descending density order
end

frags_original_order = extract_json('arabidopsis_datasets/'+ARGV[0].to_s+'/frag_ids_original_order.json')
frags_reverse_order = frags_original_order.reverse
frags_by_density = extract_json('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/frags_by_density.json')

id_pos_hash = extract_json('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/pos_hash.json') #remember the fragments are in a random order here - the order they were in in the fasta
fasta_lengths = extract_json('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/fasta_lengths.json') #also in the order from the fasta file

id_density_hash = extract_json('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/id_density_hash.json')

hs = score(frags_reverse_order, frags_original_order)
puts 'Highest possible score: ' + hs.to_s
ds = score(frags_original_order, frags_by_density)
puts 'Density order Score: ' + ds.to_s
random_scores = random_score(frags_original_order, frags_by_density)
rs = random_scores.inject(:+)/random_scores.length
puts 'Random Score: ' + rs.to_s
eos = score(frags_original_order, (even_odd(frags_by_density, 'even')))
puts 'Even Odd Method Score: ' + eos.to_s
oes = score(frags_original_order, (even_odd(frags_by_density, 'odd')))
puts 'Odd Even Method Score: ' + oes.to_s
lrs = score(frags_original_order, (left_right(id_pos_hash, fasta_lengths).flatten))
puts 'Left Right Method Score: ' + lrs.to_s
lrds = score(frags_original_order, (lr_d(frags_by_density, fasta_lengths, id_pos_hash).flatten))
puts 'Left Right Density Method Score: ' + lrds.to_s

scores = [hs, ds, rs, eos, oes, lrs, lrds]

# Rearranged orders of densities: for scatter plots
#-----------------------------------------------------------

Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/"+ARGV[0].to_s+"/re_files/id_n_density_txts"))

#Original order
scatta_txts(re_order_densities(id_density_hash, frags_original_order), 'd_o')

#Density C1
scatta_txts(re_order_densities(id_density_hash, frags_by_density), 'd_c1')

#Random C2 (example)
id_c2 = frags_by_density.shuffle
scatta_txts(re_order_densities(id_density_hash, id_c2), 'd_c2')

#Even Odd M1 a/b
id_m1a = even_odd(frags_by_density, 'even')
scatta_txts(re_order_densities(id_density_hash, id_m1a), 'd_m1a')

id_m1b = even_odd(frags_by_density, 'odd')
scatta_txts(re_order_densities(id_density_hash, id_m1b), 'd_m1b')

#Left Right M2 a: l, r, both
lnr = left_right(id_pos_hash, fasta_lengths)

id_m2al = lnr[0]
scatta_txts(re_order_densities(id_density_hash, id_m2al), 'd_m2al')

id_m2ar = lnr[1]
scatta_txts(re_order_densities(id_density_hash, id_m2ar), 'd_m2ar')

id_m2a = lnr.flatten
scatta_txts(re_order_densities(id_density_hash, id_m2a), 'd_m2a')

# Left Right Density M2 b: l, r, both
lnrd = lr_d(frags_by_density, fasta_lengths, id_pos_hash)

id_m2bl = lnrd[0]
scatta_txts(re_order_densities(id_density_hash, id_m2bl), 'd_m2bl')

id_m2br = lnrd[1]
scatta_txts(re_order_densities(id_density_hash, id_m2br), 'd_m2br')

id_m2b = lnrd.flatten
scatta_txts(re_order_densities(id_density_hash, id_m2b), 'd_m2b')


write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/random_scores.txt', random_scores)
write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/re_files/scores.txt', scores)