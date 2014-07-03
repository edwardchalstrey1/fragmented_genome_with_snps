#encoding: utf-8

require_relative '../lib/fitness_score'
require 'test/unit'

class TestFitnessScore < Test::Unit::TestCase

	def setup
		@genome_length = 30
		@div = 6
		@hm = [5, 10,13,14 ,18,19, 20, 26]
		@ht = [1,2,3,4, 6,9, 11,12, 15,16,17, 21]

		@hm_count = FitnessScore.count(@hm, @div, @genome_length)
		@ht_count = FitnessScore.count(@ht, @div, @genome_length)

		@count_ratios = FitnessScore.ratio(@hm, @ht, @div, @genome_length)
	end

	def test_count
		assert_equal([1,1,2,3,0,1], @hm_count)
		assert_equal([4,2,3,2,1,0], @ht_count)
	end

	def test_ratio
		assert_equal([2.0/5.0, 2.0/3.0, 0.75, 4.0/3.0, 0.5, 2.0], @count_ratios)
	end

	def test_score
		expected = [1,2,3,4,5] # example of an expected count ratio
		permutation = [2,1,3,4,5] # example of a count ratio from a permutation
		assert_equal(0.9, ('%.1f' % FitnessScore.score(expected, permutation)).to_f)
	end

	def test_distance_score
		expected = 5+3+1+4+1+1+6
		assert_equal(expected, FitnessScore.distance_score(@hm))
	end


end
