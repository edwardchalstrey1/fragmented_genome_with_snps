#encoding: utf-8
require_relative 'lib/write_it'
require_relative 'lib/locate_mutation'

dataset = ARGV[0]
run = ARGV[1]
gen = ARGV[2]

Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}")) do
	
	hm = WriteIt.file_to_ints_array("hm_snps.txt")
	# ht = WriteIt.file_to_ints_array("ht_snps.txt")
	puts LocateMutation.find_peak(hm)

	perm_hm = WriteIt.file_to_ints_array("p_run#{run}/Gen#{gen}_lists/gen_#{gen}_hm.txt")

	puts LocateMutation.find_peak(perm_hm)

end