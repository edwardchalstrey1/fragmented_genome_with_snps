#!/usr/env ruby

require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'pp'

vcfs = []
File.open('snps.vcf', "r").each do |line| #get array of vcf lines, you can call a method on one line
	next if line =~ /^#/
	v = Bio::DB::Vcf.new(line)
	vcfs << v.chrom
end
snps_hash = Hash.new(0)
vcfs.each {|v| snps_hash[v] +=1 } #we have the number of snps on each frag, by counting the repeats of each frag in the vcf
#the frag_id(.chrom) is the key, the number of snps for that frag is the value

fasta = [] #we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
Bio::FastaFormat.open('frags.fasta').each do |i| #get array of fasta format frags
	fasta << i
end
fasta_ids = []
fasta.each do |i|
	fasta_ids << i.entry_id
end
snps_per_frag = []
fasta_ids.each do |id|
	snps_per_frag << snps_hash[id] #gives 0 for keys that don't exist = good
end
#now we have an array with the number of snps per frag in the same order as the fasta array, which we can get lengths from to calculate density
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

### METHOD 1: frags_by_density SCORE: 398,988	density order
### METHOD 2: rearranged       SCORE: 585,252	even/odd method
### METHOD 3: rearranged2      SCORE: 403,162	somehow worse

#rearranged = frags_by_density.values_at(* frags_by_density.each_index.select {|i| i.even?}) #the even numbers first half
#rl = frags_by_density.values_at(* frags_by_density.each_index.select {|i| i.odd?}) #the odd numbers reversed second half
#rl.reverse_each do |i|
#	rearranged << i        #THIS METHOD OF REARRANGEMENT IS SHIT, mostly because we know that a lot of the 0 density frags are in order
#end
#puts rearranged.length

rearranged2 = ids_zero_snps
lrs = snp_ids_density_order.values_at(* snp_ids_density_order.each_index.select {|i| i.even?})
rls = snp_ids_density_order.values_at(* snp_ids_density_order.each_index.select {|i| i.odd?})
rearranged2 << rls
lrs.reverse_each do |i|
	rearranged2 << i
end
rearranged2 << ids_zero_snps_end
rearranged2 = rearranged2.flatten
puts rearranged2.length

###ALL THE CODE BELOW USED TO DETERMINE SIMILARITY BETWEEN ORIGINAL AND REARRANGED FRAGS
###THIS WORKS BECAUSE THE FRAGS IN THE FASTA ARE IN ORDER! #will need to make a version that doesn't rely on this eventually
position_each_frag_id_in_d = fasta_ids.map{|x| rearranged2.index(x)} #works out the index of fasta_id values in frags_by_density
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


