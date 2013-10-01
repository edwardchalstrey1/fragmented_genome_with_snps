#!/usr/env ruby

require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'pp'

vcfs = []
File.open('snps.vcf', "r").each do |line| #get array of vcf lines, you can call a method on one line
	next if line =~ /^#/
	v = Bio::DB::Vcf.new(line)
	vcfs << v
end
snps_hash = Hash.new(0)
vcfs.each {|v| snps_hash[v.chrom] +=1 } #we have the number of snps on each frag, by counting the repeats of each frag in the vcf
#the frag_id(.chrom) is the key, the number of snps for that frag is the value

fasta = [] #we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
Bio::FastaFormat.open('frags.fasta').each do |i| #get array of fasta format frags
	fasta << i
end
fasta_ids = []
fasta.each do |i|
	fasta_ids << i.entry_id
end
snps_per_frag = [] #use the id to identify the number of snps for that frag using the keys of snps_hash
fasta_ids.each do |id|
	snps_per_frag << snps_hash[id] #gives 0 for keys that don't exist = good, because the frags with 0 density would otherwise be missing
end
#now we have an array with the number of snps per frag in the same order as the fasta array, which we can get lengths from to calculate density
vcfs_chrom = []
vcfs.each do |v|
	vcfs_chrom << v.chrom
end
vcfs_pos = []
vcfs.each do |v|
	vcfs_pos << v.pos
end

pos = [] #get the positions for each frag, in an array of arrays
n = 0
fasta_ids.each do |i|
	each_fr_pos = []
	while i == vcfs_chrom[n]
		each_fr_pos << vcfs_pos[n]
		n+=1
	end
	pos << each_fr_pos
	n = pos.flatten.length
end
fasta_lengths = []
fasta.each do |i|
	fasta_lengths << i.length
end
#make hash where key is frag_id and values are positions for each frag, positions array should be converted to string so they stick together
#otherwise we have problems when flattening the array to make the hash below
pos_strings = []
pos.each do |s|
	pos_strings << "#{s.join(",")}" #each position separated by a comma
end
pos_hash = Hash[*fasta_ids.zip(pos_strings).flatten] # fasta ids as keys and the snp position strings as values

left = [] # frags with snps's mostly to the right, so should be on left of distribution
right = [] # frags with snps's mostly to the left, so should be on right of distribution
d = 0
fasta_ids.each do |x| # for each of the fragments
	temp = []
	(pos_hash[x].split(",").map { |s| s.to_i }).each do |i| # adding each of the positions associated with that frag/key to temporaray array
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
	end
	d+=1 
end

densities = []
ids_w_snps = []
ids_zero_snps = []
ids_zero_snps_end = []
x = 0
fasta.each do |frag| #find a way to split the first half of the zero snps from the second half (ends of normal dist)
	if snps_per_frag[x] == 0 && ids_w_snps.length == 0
		ids_zero_snps << frag.entry_id.to_s
		x+=1
	elsif snps_per_frag[x] != 0
		densities << (snps_per_frag[x].to_f / frag.length.to_f)*1000 #this gives snps/Kb as density units
		ids_w_snps << frag.entry_id.to_s
		x+=1
	elsif snps_per_frag[x] == 0 && ids_w_snps.length != 0
		ids_zero_snps_end << frag.entry_id.to_s
		x+=1
	end
end
a = densities.zip(ids_w_snps) # array of density and id values, for the frags w snps
density_order = a.sort 
snp_ids_density_order = density_order.flatten.values_at(* density_order.flatten.each_index.select {|i| i.odd?})
# selecting the odd values from the flattened density order array produces a list of frag ids ranked by density
frags_by_density = [] 
frags_by_density << ids_zero_snps #all of these frags have 0 snp density
frags_by_density << ids_zero_snps_end
frags_by_density << snp_ids_density_order
frags_by_density = frags_by_density.flatten #all the frag ids now ranked by density with all the zero's at the start

### METHOD 1: frags_by_density   SCORE: 385,572	density order
### METHOD 2: rearranged         SCORE: 520,762	even/odd method
### METHOD 3: rearranged2        SCORE: 166,972	somehow worse - better, was using shuffled fasta before
### METHOD 4: left_right_method  SCORE: 155,496	SHOULD REARRANGE THE LEFT AND RIGHT BY DENSITY!!
### METHOD 5: left_right_method2 SCORE: 168,178 surprisingly worse

#make array of order 0,left,right,0end
left_right_method = []
left_right_method << ids_zero_snps
left_right_method << left
left_right_method << right
left_right_method << ids_zero_snps_end
left_right_method = left_right_method.flatten
#puts left_right_method.length

aa = ids_w_snps.zip(densities)
aa_hash = Hash[*aa.flatten]
left_densities = []
left.each do |i|
	left_densities << aa_hash[i] #finding the densities associated with the frags in left
end
right_densities = []
right.each do |i|
	right_densities << aa_hash[i] #finding the densities associated with the frags in right
end
ll = left_densities.zip(left) # getting the 'left' and 'right' frags into density order
left_density_order = ll.sort.flatten.values_at(* ll.sort.flatten.each_index.select {|i| i.odd?})
rr = right_densities.zip(right)
right_density_order = rr.sort.flatten.values_at(* rr.sort.flatten.each_index.select {|i| i.odd?})
right_density_order.reverse! # the order of the 'right' frags needs to be density descending

left_right_method2 = []
left_right_method2 << ids_zero_snps
left_right_method2 << left_density_order
left_right_method2 << right_density_order
left_right_method2 << ids_zero_snps_end
left_right_method2 = left_right_method2.flatten
#puts left_right_method2.length

###ALL THE CODE BELOW USED TO DETERMINE SIMILARITY BETWEEN ORIGINAL AND REARRANGED FRAGS
###THIS WORKS BECAUSE THE FRAGS IN THE FASTA ARE IN ORDER! #will need to make a version that doesn't rely on this eventually
position_each_frag_id_in_d = fasta_ids.map{|x| left_right_method2.index(x)} #works out the index of fasta_id values in frags_by_density
index_values = Array(0..(fasta_ids.length - 1)) # index values that fasta_ids originally at
both = []
both << position_each_frag_id_in_d
both << index_values
difference = both.transpose.map {|x| x.reduce(:-)} #these differences will obviously be a lot, a first experiment
# taking away old position from new position, to find the distance that the frag has moved when re-ordered by density
# when doing the method above again for the reordered frags_by_density, all values will be 0, if order the same as frag_ids
difference_abs = []
difference.each do |i|
	difference_abs << i.abs
end
largest_diff = difference_abs.each_with_index.max[0]
differences_proportions = []
difference_abs.each do |i|
	p = i.to_f/largest_diff.to_f
	differences_proportions << p #this proportional score could be useful in linear programming
end
#puts differences_proportions #the closer to 1, the more wrong. proportion relative to the one furthest from original position
score = difference_abs.inject(:+) #high score = bad, score of 0 means the fragments in the right order
puts score
#puts frags_by_density


