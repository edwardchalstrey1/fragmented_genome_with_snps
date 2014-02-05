require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'
require_relative '~/fragmented_genome_with_snps/lib/reform_ratio'
require 'test/unit'


# UNIT TESTING THE ReformRatio CLASS

class TestReform < Test::Unit::TestCase
	def test_1
		a = ["a", "b", "c"]
		b = ["a", "b", "c"]
		c = ["c", "b", "a"]
		assert_equal(0, ReformRatio::rearrangement_score(a,b))
		assert_equal(4, ReformRatio::rearrangement_score(c,b))
	end
end