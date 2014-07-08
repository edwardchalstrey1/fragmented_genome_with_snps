#encoding: utf-8
require_relative '../lib/example_perms'
require_relative '../lib/reform_ratio'
require 'test/unit'

class FakeFasta
	attr_accessor :entry_id, :length
	def initialize
		@entry_id = 'FragX'
		@length = 10
	end
end

class TestExample < Test::Unit::TestCase

	def fasta
		frag1, frag2, frag3, frag4 = FakeFasta.new, FakeFasta.new, FakeFasta.new, FakeFasta.new
		frag1.entry_id, frag2.entry_id, frag3.entry_id, frag4.entry_id = 'frag1', 'frag2', 'frag3', 'frag4'
		return [frag1, frag2, frag3, frag4], [frag4, frag3, frag2, frag1]
	end

	def test_fasta_p_id
		fasta1, fasta2 = fasta
		fasta2_ids = ReformRatio.fasta_id_n_lengths(fasta2)[0]
		fasta1_reordered = ExamplePerms.fasta_p_id(fasta1, fasta2_ids)
		assert_equal(fasta2, fasta1_reordered)
	end

	def test_get_perms
		all_perms = ExamplePerms.get_perms(2, 0, 1, 'small_dataset2', 'table_test')
		assert_equal(2, all_perms.length)
		assert_kind_of(Array, all_perms[0], 'population not array')
		assert_equal(10, all_perms[0].length, 'wrong number of permutations')
		assert_equal(54, all_perms[0][0].length, 'permutation wrong length') # frags plus fitness
		assert_kind_of(String, all_perms[0][0][0])
	end
end