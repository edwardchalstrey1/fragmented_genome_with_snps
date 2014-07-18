#encoding: utf-8
require_relative 'lib/write_it'
require_relative 'lib/locate_mutation'
require_relative 'lib/snp_dist'
require_relative 'lib/reform_ratio'
require_relative 'lib/fitness_score'

dataset = ARGV[0]
run = ARGV[1]
gen = ARGV[2]
n = 1048576*4

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")
div = 10000

Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}")) do
	
	hm = WriteIt.file_to_ints_array("hm_snps.txt") # Get SNP distributions
	ht = WriteIt.file_to_ints_array("ht_snps.txt")

	ratios = FitnessScore::ratio(hm, ht, div, genome_length) # Calculate homozygous/heterozygous ratio and make approximate distribution
	hyp = SNPdist.hyp_snps(ratios, genome_length)

	peak =  LocateMutation.find_peak(hyp, n) # Find the peak in the approximated (hypothetical SNP) distribution
	causal = LocateMutation.closest_snp(peak, hm)

	perm_hm = WriteIt.file_to_ints_array("p_run#{run}/Gen#{gen}_lists/gen_#{gen}_hm.txt")
	perm_ht = WriteIt.file_to_ints_array("p_run#{run}/Gen#{gen}_lists/gen_#{gen}_ht.txt")

	perm_ratios = FitnessScore::ratio(perm_hm, perm_ht, div, genome_length) # Calculate homozygous/heterozygous ratio and make approximate distribution
	perm_hyp = SNPdist.hyp_snps(perm_ratios, genome_length)

	perm_peak = LocateMutation.find_peak(perm_hyp, n)
	candidate = LocateMutation.closest_snp(perm_peak, perm_hm)

	puts causal
	puts candidate

end