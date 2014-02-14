#encoding: utf-8
require_relative 'lib/rearrangement_score'
require 'test/unit'

class TestRearrangementScore < Test::Unit::TestCase

	def test_rearrangement_score
		a = %w(a b c)
		c = a.reverse

		assert_equal(0, RearrangementScore::rearrangement_score(a,a))
		assert_equal(4, RearrangementScore::rearrangement_score(c,a))

		assert_equal('0/4', RearrangementScore::score(a,a))
	end
end