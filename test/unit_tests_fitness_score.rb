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

		@ratios = [2.0/5.0, 2.0/3.0, 0.75, 4.0/3.0, 0.5, 2.0]

		@snps1 = [1,5,6,8,12]
		@snps2 = [1,11,21,31,41,51]
	end

	def test_count
		assert_equal([1,1,2,3,0,1], @hm_count)
		assert_equal([4,2,3,2,1,0], @ht_count)
	end

	def test_ratio
		assert_equal(@ratios, FitnessScore.ratio(@hm, @ht, @div, @genome_length))
	end

	def test_count_ratio
		assert_equal(1.0, FitnessScore.count_ratio(@hm, @ht, @div, @genome_length, @ratios).round(2))
	end

	def test_snp_distance
		score = FitnessScore.snp_distance(@snps1)
		assert_equal(11, score)
	end

	def test_max_density
		score = FitnessScore.max_density(@snps2)
		assert_equal(0.01654088, ('%.8f' % score).to_f)
	end

	def test_max_ratio
		score = FitnessScore.max_ratio(@snps1, @snps2)
		assert_equal(15.37267, ('%.5f' % score).to_f)
	end

	def test_max_hyp
		score = FitnessScore.max_hyp(@snps1, @snps2, 6, 51)
		assert_kind_of(Float, score)
	end
end
