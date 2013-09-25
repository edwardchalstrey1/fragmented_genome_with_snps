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
vcfs.each {|v| snps_hash[v] +=1 } #we have the number of snps on each frag

fasta = [] #we have the lengths of each fasta, but they aren't in the same order as in the vcf/hash
Bio::FastaFormat.open('frags.fasta').each do |i| #get array of fasta format frags
	fasta << i
end
fasta_ids = []
fasta.each do |i|
	fasta_ids << i.entry_id
end
snps_per_frag = []
fasta.each do |frag|
	snps_per_frag << snps_hash[frag.entry_id]
end
#now we have an array with the number of snps per frag in the same order as the fasta array, which we can get lengths from to calculate density
densities = []
x = 0
fasta.each do |frag|
	if snps_per_frag[x] == 0
		x +=1
	else 
		densities << (snps_per_frag[x].to_f / frag.length.to_f)*1000 #this gives snps/Kb as density units
		x += 1
	end
end
a = densities.zip(fasta_ids).flatten.compact
density_frag_hash = Hash[a.each_slice(2).to_a]
sorted_density_frag_hash = density_frag_hash.sort
puts sorted_density_frag_hash
