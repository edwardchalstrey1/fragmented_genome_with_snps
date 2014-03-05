#encoding: utf-8
require_relative 'lib/rearrangement_score'
require 'test/unit'

class TestRearrangementScore < Test::Unit::TestCase

	A = %w(a b c)
	C = A.reverse

	def test_dev_dist
		assert_equal(0, RearrangementScore::dev_dist(A,A))
		assert_equal(1, RearrangementScore::dev_dist(A,C), 'ac dev dist wrong')
		assert_equal(1, RearrangementScore::dev_dist(%w(a b), %w(b a)))
	end

	def test_mod_ham_dist
		assert_equal(0, RearrangementScore::mod_ham_dist(A,A))
		assert_equal(1, RearrangementScore::mod_ham_dist(A,C))
		assert_equal(5.0/6.0, RearrangementScore::mod_ham_dist(%w(a b c d), %w(d c b a)))
		assert_equal(1, RearrangementScore::mod_ham_dist(%w(a b c d), %w(d b c a)))
	end

	def test_gen_ham_dist
		assert_equal(0, RearrangementScore::gen_ham_dist(A,A))
		assert_equal(1, RearrangementScore::gen_ham_dist(A, %w(b c a)))
	end

	def test_sq_dev_dist
		assert_equal(0, RearrangementScore::sq_dev_dist(A,A))
		assert_equal(1, RearrangementScore::sq_dev_dist(A,C))
	end

	def test_r_dist
		assert_equal(0, RearrangementScore::r_dist(A,A))
		assert_equal(1, RearrangementScore::r_dist(A,C))
	end

	def test_lcs
		assert_equal(0, RearrangementScore::lcs(A,A))
		assert_equal(1, RearrangementScore::lcs(A,C))
	end

	def test_kendalls_tau
		assert_equal(0, RearrangementScore::kendalls_tau(A,A))
		assert_equal(1, RearrangementScore::kendalls_tau(A,C))
	end
end