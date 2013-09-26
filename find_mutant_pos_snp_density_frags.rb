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
a = densities.zip(ids_w_snps)
density_order = a.sort
snp_ids_density_order = density_order.flatten.values_at(* density_order.flatten.each_index.select {|i| i.odd?})
frags_by_density = []
frags_by_density << ids_zero_snps #all of these frags have 0 snp density
frags_by_density << snp_ids_density_order
frags_by_density = frags_by_density.flatten

x = 0 #the value of x will represent the position (index) in the density array
position_each_frag_id_in_d = [] #want to get positions of the values in frag_ids in frags_by_density
iteration = []
fasta_ids.each do |i|
	if frags_by_density[x] == i
		position_each_frag_id_in_d << x #if the value at position x matches the value at i, add it to the new array
		iteration << i
	else
		until frags_by_density[x] == i #otherwise increment x until they do match, and add the position
			x +=1
		end
		position_each_frag_id_in_d << x
		iteration << i
	end
	x = iteration.length # x should be incremented, however I cannot simply do: x += 1, as x may have been incremented by the until loop
	puts x
end
puts position_each_frag_id_in_d

