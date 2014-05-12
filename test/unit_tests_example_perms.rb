#encoding: utf-8
require_relative '../lib/score_plots/example_perms'
require 'bio'
require 'test/unit'

class TestExample < Test::Unit::TestCase

	Fasta = ReformRatio::fasta_array("arabidopsis_datasets/small_dataset2/frags.fasta")

	def test_fasta_p
		permutation = ExamplePerms::fasta_p(Fasta, "arabidopsis_datasets/small_dataset2/run13/Gen12/best_permutation.txt")
		assert_kind_of(Array, permutation)
		assert_kind_of(Float, permutation[0])
		assert_in_delta(0.5, permutation[0], 0.5)
		assert_equal(53, permutation[1..-1].length)
		assert_kind_of(Bio::FastaFormat, permutation[1])
	end
end