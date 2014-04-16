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

