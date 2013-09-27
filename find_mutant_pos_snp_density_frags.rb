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
x = 0
fasta.each do |frag|
	if snps_per_frag[x] == 0
		ids_zero_snps << frag.entry_id.to_s
		x +=1
	else 
		densities << (snps_per_frag[x].to_f / frag.length.to_f)*1000 #this gives snps/Kb as density units
		ids_w_snps << frag.entry_id.to_s
		x += 1
	end
end
a = densities.zip(ids_w_snps) # array of density and id values, for the frags w snps
density_order = a.sort 
snp_ids_density_order = density_order.flatten.values_at(* density_order.flatten.each_index.select {|i| i.odd?})
# selecting the odd values from the flattened density order array produces a list of frag ids ranked by density
frags_by_density = [] 
frags_by_density << ids_zero_snps #all of these frags have 0 snp density
frags_by_density << snp_ids_density_order
frags_by_density = frags_by_density.flatten #all the frag ids now ranked by density with all the zero's at the start

position_each_frag_id_in_d = fasta_ids.map{|x| frags_by_density.index(x)} #works out the index of fasta_id values in frags_by_density
index_values = Array(0..(fasta_ids.length - 1)) # index values that fasta_ids originally at
both = []
both << position_each_frag_id_in_d
both << index_values
difference = both.transpose.map {|x| x.reduce(:-)} #these differences will obviously be a lot, a first experiment
# taking away old position from new position, to find the distance that the frag has moved when re-ordered by density
# when doing the method above again for the reordered frags_by_density, all values will be 0, if order the same as frag_ids
puts difference #DISTANCE FROM ORIGINAL POSITION - positive no. means moved forward in array, negative back




