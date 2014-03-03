#encoding: utf-8
require_relative 'lib/rearrangement_score'
require 'test/unit'

class TestRearrangementScore < Test::Unit::TestCase

	A = %w(a b c)
	C = A.reverse

	def test_rearrangement_score
		assert_equal(0, RearrangementScore::rearrangement_score(A,A))
		assert_equal(4, RearrangementScore::rearrangement_score(A,C))

		assert_equal('0/4', RearrangementScore::score(A,A))
	end

	def test_metric_2
		assert_equal(0, RearrangementScore::metric_2(A,A, 'a'))
		assert_equal(4, RearrangementScore::metric_2(A,C, 'a'))

		assert_equal(2+(2.0/3.0), RearrangementScore::metric_2(A,C, 'b'))
		assert_equal(2+(2.0/3.0), RearrangementScore::metric_2(A,C, 'c'))
	end

	def test_metric_1_2_av
		assert_equal(4, RearrangementScore::metric_1_2_av(A,C))
	end

end