#encoding: utf-8
require_relative 'lib/write_it'
require_relative 'lib/locate_mutation'
require_relative 'lib/snp_dist'
require_relative 'lib/reform_ratio'
require_relative 'lib/fitness_score'

dataset = ARGV[0]
run = ARGV[1]
gen = ARGV[2]

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")

Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}")) do
	
	hm = WriteIt.file_to_ints_array("hm_snps.txt")
	ht = WriteIt.file_to_ints_array("ht_snps.txt")
	comparable_ratio = FitnessScore::ratio(hm, ht, 100, genome_length)
	hyp = SNPdist.hyp_snps(comparable_ratio, genome_length)

	puts LocateMutation.find_peak(hm, 1024)#1048576*2)

	# perm_hm = WriteIt.file_to_ints_array("p_run#{run}/Gen#{gen}_lists/gen_#{gen}_hm.txt")

	# puts LocateMutation.find_peak(perm_hm, n)

end