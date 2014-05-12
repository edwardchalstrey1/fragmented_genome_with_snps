#encoding: utf-8
require_relative 'lib/score_plots/example_perms'
require_relative 'lib/snp_dist'
require_relative 'lib/reform_ratio'

=begin
	
ARGV's: 0 = dataset
		1 = run
		2 = generation (number)
		3 = number of divisions of the genome where the ratio of homozygous to heterozygous SNPs is calculated
=end

div = ARGV[3].to_f

original_fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")
permutation_txt = "arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{ARGV[2]}/best_permutation.txt"
fasta = ExamplePerms::fasta_p(original_fasta, permutation_txt)[1..-1]

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf")

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")
snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) # array of no. of snps per frag in same order as fasta
pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) # get snp positions for each frag in array of arrays
actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
het_snps, hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])

fratio_breaks_perm = SNPdist::fratio(hom_snps, het_snps, div, genome_length) # frequency ratio array
perm_ratio = SNPdist::hyp_snps(fratio_breaks_perm, div, genome_length) # hypothetical snp positions array

SNPdist::plot_hyp(perm_ratio, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", ARGV[2], genome_length)
