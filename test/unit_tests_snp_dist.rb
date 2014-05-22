#encoding: utf-8
require_relative '../lib/snp_dist'
require_relative '../lib/model_genome'
require 'test/unit'

class TestSNPdist < Test::Unit::TestCase

	def setup
		@div = 100.0; @genome_length = 2000.0
		hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
		@fratio, @breaks = SNPdist::fratio(hm, ht, @div, @genome_length)
		@comparable_ratio = SNPdist::hyp_snps([@fratio, @breaks], @div, @genome_length)
	end

	def test_fratio
		assert_equal(@div, @fratio.length)
		assert_equal(@div+1, @breaks.length)
		assert_kind_of(Integer, @fratio[0])
		assert_kind_of(Float, @breaks[1])
	end

	def test_hyp_snps
		assert_kind_of(Array, @comparable_ratio)
		assert_kind_of(Float, @comparable_ratio[0])
	end

	def test_qq_cor
		hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
		fratio_breaks = SNPdist::fratio(hm, ht, @div, @genome_length)
		ratio = SNPdist::hyp_snps(fratio_breaks, @div, @genome_length)
		cor = SNPdist::qq_cor(@comparable_ratio, ratio)
		assert_kind_of(Float, cor)
		assert(0 <= cor && cor <= 1)
	end

end