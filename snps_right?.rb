#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require 'pp'

dataset = ARGV[0]
location = 'fragmented_genome_with_snps/arabidopsis_datasets'
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"

snp_data = ReformRatio::get_snp_data(vcf_file)
fasta = ReformRatio::fasta_array(fasta_file)

snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) # array of no. of snps per frag in same order as fasta
pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) # get snp positions for each frag in array of arrays
actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
het_snps, hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])

hm = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/hm_snps.txt") # we can use the SNPs from the model genome to make example ratio
ht = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/ht_snps.txt")

puts het_snps.sort == ht.sort
puts hom_snps.sort == hm.sort

# pp het_snps.sort
# pp ht.sort

het_snps.each do |snp|
	unless ht.include?(snp)
		puts "ht dunt hav: #{snp}"
	end
end
ht.sort.each do |snp|
	unless het_snps.include?(snp)
		puts "het dunt hav: #{snp}"
	end
end

puts het_snps.length
puts ht.length