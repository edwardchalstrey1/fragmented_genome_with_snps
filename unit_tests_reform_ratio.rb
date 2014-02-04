require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'
require_relative 'lib/reform_ratio'
require 'test/unit'


# UNIT TESTING THE ReformRatio CLASS

class FakeFasta
	attr_accessor :entry_id, :length
	def initialize
		@entry_id = "FragX"
		@length = 10
	end
end

class TestReform < Test::Unit::TestCase
	def test_rearrangement_score
		a = ["a", "b", "c"]
		b = ["a", "b", "c"]
		c = ["c", "b", "a"]
		assert_equal(0, ReformRatio::rearrangement_score(a,b))
		assert_equal(4, ReformRatio::rearrangement_score(c,b))
	end
	def test_total_pos
		pos = [[1,2,3],[4,5,6],[7,8,9]]
		lengths = [10,10,10]
		expected = [1,2,3,13,14,15,26,27,28]
		assert_equal(expected, ReformRatio::total_pos(pos, lengths))
	end
	def test_prime?
		assert_equal(true, ReformRatio::prime?(7))
		assert_equal(false, ReformRatio::prime?(4))
		assert_equal(true, ReformRatio::prime?(2))
		assert_equal(true, ReformRatio::prime?(971))
		assert_equal(false, ReformRatio::prime?(972))
	end
	def test_division
		a = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
		b = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
		ax = ReformRatio::division(a)
		bx = ReformRatio::division(a)
		assert(ax==2||ax==1||ax==4||ax==5||ax=10||ax==20) 
		assert(bx==18||bx==9||bx=6||bx==3||bx==2||bx==1)
	end
	def test_recombination
		parent1 = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j","k", "l", "m", "n", "o", "p", "q", "r", "s", "t"] #20
		parent2 = parent1.reverse
		child = ReformRatio::recombine(parent1, parent2)

		parent3 = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j","k", "l", "m", "n", "o", "p", "q", "r", "s"] #19
		parent4 = parent3.reverse
		child2 = ReformRatio::recombine(parent3, parent4)
		
		assert(child.uniq == child, "Child of p1/2 not unique")
		assert(child != parent1, "Child same as parent1")
		assert(child != parent2, "Child same as parent2")

		assert(child2.uniq == child2, "Child of p3/4 not unique")
		assert(child2 != parent3, "Child same as parent3")
		assert(child2 != parent4, "Child same as parent4")
	end
	def test_fasta_id_n_lengths
		fragment = FakeFasta.new
		fasta = [fragment, fragment, fragment]
		ids_n_lengths = ReformRatio::fasta_id_n_lengths(fasta)
		ids = ids_n_lengths[0]
		lengths = ids_n_lengths[1]
		assert_equal(['FragX', 'FragX', 'FragX'], ids)
		assert_equal([10,10,10], lengths)
	end
	def test_snps_per_fasta_frag
		frag1, frag2 = FakeFasta.new, FakeFasta.new
		frag1.entry_id, frag2.entry_id = 'frag1', 'frag2'
		fasta = [frag1, frag2]
		h = {'frag2'=>3,'frag1'=>2}
		assert_equal([2,3], ReformRatio::snps_per_fasta_frag(h, fasta))
	end
end

