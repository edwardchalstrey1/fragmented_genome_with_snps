require 'bio-samtools'
require 'bio'
require 'pp'
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'

perm_ids = []
IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{ARGV[2]}_best_permutation.txt") { |line| perm_ids << line.gsub(/\n/,'') }
perm_ids = perm_ids[2..-1]

fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")

ids = ReformRatio::fasta_id_n_lengths(fasta)[0]

fasta_perm = []
perm_ids.each do |id|
	fasta.each do |frag|
		if id == frag.entry_id
			fasta_perm << frag
		end
	end
end

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf")
snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta_perm) #array of no. of snps per frag in same order as fasta
pos_n_info = ReformRatio::get_positions(fasta_perm, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) #get snp positions for each frag in array of arrays
actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta_perm)[1])
het_hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])

WriteIt::write_txt("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/gen#{ARGV[2]}_best_permutation_het_snps", het_hom_snps[0]) #het
WriteIt::write_txt("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/gen#{ARGV[2]}_best_permutation_hom_snps", het_hom_snps[1]) #hom