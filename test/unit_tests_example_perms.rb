#encoding: utf-8
require_relative '../lib/score_plots/example_perms'
require_relative '../lib/model_genome'
require_relative '../lib/snp_dist'
require 'bio'
require 'test/unit'

class TestExample < Test::Unit::TestCase

	def setup
		@div = 100.0; @genome_length = 2000.0
		hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
		fratio, breaks = SNPdist::fratio(hm, ht, @div, @genome_length)
		@hyp = SNPdist::hyp_snps([fratio, breaks], @div, @genome_length)

		@fasta = ReformRatio::fasta_array("arabidopsis_datasets/small_dataset2/frags.fasta")
		@snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/small_dataset2/snps.vcf")
	end

	def test_fasta_p
		permutation = ExamplePerms::fasta_p(@fasta, "arabidopsis_datasets/small_dataset2/run13/Gen12/best_permutation.txt")
		assert_kind_of(Array, permutation)
		assert_kind_of(Float, permutation[0])
		assert_in_delta(0.5, permutation[0], 0.5)
		assert_equal(53, permutation[1..-1].length)
		assert_kind_of(Bio::FastaFormat, permutation[1])
	end

	def test_get_perms
		perms = ExamplePerms::get_perms(@fasta, 2, @snp_data, @hyp, @div, @genome_length)
		assert_equal(3, perms.length)
		perms.each do |pop|
			assert_kind_of(Array, pop)
			assert_kind_of(Array, pop[0])
			assert_kind_of(Float, pop[0][0])
			assert_kind_of(String, pop[0][1], 'should be a frag_id string')
			assert_kind_of(String, pop[0][-1], 'should be frag_id')
			assert_equal(54, pop[0].length)
		end
	end
end