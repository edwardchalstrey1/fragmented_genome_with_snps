#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'
require_relative 'lib/model_genome'
require_relative 'lib/snp_dist'

if ARGV[0] == nil || ARGV[1] == nil
	puts 'Wrong number of command line arguments! Specify a dataset and test-set'
else

	vcf_file = "arabidopsis_datasets/#{ARGV[0]}/snps.vcf"
	fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags_shuffled.fasta"

	hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
	fratio_breaks = SNPdist::fratio(hm, ht, 10000) # frequency ratio array
	comparable_ratio = SNPdist::hyp_snps(fratio_breaks, 10000) # hypothetical snp positions array

	GATOC::evolve(fasta_file, vcf_file, :pop_size => 20, :select_num => 10, :mut_num => 2,
	 :save => 1, :ran => 1, :loc => 'fragmented_genome_with_snps/arabidopsis_datasets', :comparable_ratio => comparable_ratio)
end
