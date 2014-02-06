#encoding: utf-8
require 'rubygems'
require 'rinruby'
require 'bio-samtools'
require 'bio'
require_relative 'lib/arabidopsis_c4_w_snps'
require 'test/unit'

class TestModelGenome < Test::Unit::TestCase

	def test_fasta_to_char_array
		assert_equal(%w(A T G C A T A A A A A), ModelGenome::fasta_to_char_array('test/dummy.fasta'))
	end

	def test_normal_dist
		dist = ModelGenome::normal_dist
		assert_equal(dist[0], [dist[2],dist[1]].flatten)
		assert_equal(dist[0], dist[0].uniq, "SNP pos not unique!")
	end

	def test_snp_seq
		seq = %w(A T C G N A A)
		snp_pos = [1,2,3,4,5]
		assert_equal(%w(T A G C R A A), ModelGenome::snp_seq(seq, snp_pos))
	end

	def test_get_frags
		seq = ['A']*20
		frags = ModelGenome::get_frags(seq, 2)
		assert_equal(seq.length, frags.flatten.length, 'Extra frags!')
	end
end