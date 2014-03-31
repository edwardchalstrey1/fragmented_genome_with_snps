#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'
require_relative 'lib/model_genome'
require_relative 'lib/snp_dist'
require_relative 'lib/write_it'

if ARGV[0] == nil || ARGV[1] == nil
	puts 'Wrong number of command line arguments! Specify a dataset and test-set'
else
	## Files ##
	vcf_file = "arabidopsis_datasets/#{ARGV[0]}/snps.vcf"
	fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags_shuffled.fasta"
	location = 'fragmented_genome_with_snps/arabidopsis_datasets'


	## Comparable ratio ##
	div = 100.0
	genome_length = 2000.0 # 18585056.0
	hm = WriteIt::file_to_array("#{Dir.home}/#{location}/#{ARGV[0]}/hm_snps.txt") # we can use the SNPs from the model genome to make example ratio
	ht = WriteIt::file_to_array("#{Dir.home}/#{location}/#{ARGV[0]}/ht_snps.txt")
	# hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
	fratio_breaks = SNPdist::fratio(hm, ht, div, genome_length) # frequency ratio array # 10,000 for 10K dataset
	comparable_ratio = SNPdist::hyp_snps(fratio_breaks, div, genome_length) # hypothetical snp positions array


	## Initial plots of correct ratio distribution ##
	SNPdist::plot_hyp(comparable_ratio, location, ARGV[0], 'correct', genome_length)
	SNPdist::plot_dens(hm, ht, location, ARGV[0], genome_length)


	GATOC::evolve(fasta_file, vcf_file, :pop_size => 20, :select_num => 10, :mut_num => 2,
	 :save => 1, :ran => 1, :loc => 'fragmented_genome_with_snps/arabidopsis_datasets', :comparable_ratio => comparable_ratio, :div => div, :genome_length => genome_length)
end