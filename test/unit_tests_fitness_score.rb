#encoding: utf-8

require_relative '../lib/fitness_score'
require_relative '../lib/write_it'
require 'test/unit'

class TestFitnessScore < Test::Unit::TestCase

	def setup
		@genome_length = 30
		@div = 6
		@hm = [5, 10,13,14 ,18,19, 20, 26]
		@ht = [1,2,3,4, 6,9, 11,12, 15,16,17, 21]

		@hm_count = FitnessScore.count(@hm, @div, @genome_length)
		@ht_count = FitnessScore.count(@ht, @div, @genome_length)

		@count_ratios = FitnessScore.ratio(@hm_count, @ht_count)

		@div = 100.0; @genome_length = 2000.0
		hm = WriteIt.file_to_ints_array('test/test/hm_snps.txt')
		ht = WriteIt.file_to_ints_array('test/test/ht_snps.txt')
		hom_count = FitnessScore::count(hm, @div, @genome_length)
		het_count = FitnessScore::count(ht, @div, @genome_length)
		@ratios = FitnessScore::ratio(hom_count, het_count)
	end

	def test_count
		assert_equal([1,1,2,3,0,1], @hm_count)
		assert_equal([4,2,3,2,1,0], @ht_count)
	end

	def test_ratio
		assert_equal([0.25, 0.5, 2.0/3.0, 1.5, 'NaN', 'NaN'], @count_ratios)
	end

	def test_score
		expected = [1,2,3,4,5] # example of an expected count ratio
		permutation = [2,1,3,4,5] # example of a count ratio from a permutation
		assert_equal(0.9, ('%.1f' % FitnessScore.score(expected, permutation)).to_f)

		ex2 = ['NaN',2,3,4,5]
		perm2 = [2,1,4,'NaN',5] 
		assert_equal(0.89, ('%.2f' % FitnessScore.score(ex2, perm2)).to_f)

		assert_equal(1.0, FitnessScore.score([0.25, 0.5, 2.0/3.0, 1.5, 'NaN', 'NaN'], [0.25, 0.5, 2.0/3.0, 1.5, 'NaN', 'NaN']))
		assert_equal(1.0, FitnessScore.score(@ratios, @ratios))
	end
end
