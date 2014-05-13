#encoding: utf-8
require_relative '../lib/score_plots/example_perms'
require_relative '../lib/model_genome'
require_relative '../lib/snp_dist'
require 'bio'
require 'test/unit'
require 'pp'

class TestExample < Test::Unit::TestCase

	Div = 100.0; Genome_length = 2000.0
	hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
	Fratio, Breaks = SNPdist::fratio(hm, ht, Div, Genome_length)
	Hyp = SNPdist::hyp_snps([Fratio, Breaks], Div, Genome_length)

	Fasta = ReformRatio::fasta_array("arabidopsis_datasets/small_dataset2/frags.fasta")
	Snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/small_dataset2/snps.vcf")

	def test_fasta_p
		permutation = ExamplePerms::fasta_p(Fasta, "arabidopsis_datasets/small_dataset2/run13/Gen12/best_permutation.txt")
		assert_kind_of(Array, permutation)
		assert_kind_of(Float, permutation[0])
		assert_in_delta(0.5, permutation[0], 0.5)
		assert_equal(53, permutation[1..-1].length)
		assert_kind_of(Bio::FastaFormat, permutation[1])
	end

	def test_get_perms
		perms = ExamplePerms::get_perms(Fasta, 2, Snp_data, Hyp, Div, Genome_length)
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