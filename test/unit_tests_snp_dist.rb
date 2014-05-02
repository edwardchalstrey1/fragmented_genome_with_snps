#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/snp_dist'
require '~/fragmented_genome_with_snps/lib/model_genome'
require 'test/unit'

class TestSNPdist < Test::Unit::TestCase

	Div = 100.0; Genome_length = 2000.0
	hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
	Fratio, Breaks = SNPdist::fratio(hm, ht, Div, Genome_length)
	Hyp = SNPdist::hyp_snps([Fratio, Breaks], Div, Genome_length)

	def test_fratio
		assert_equal(Div, Fratio.length)
		assert_equal(Div+1, Breaks.length)
		assert_kind_of(Integer, Fratio[0])
		assert_kind_of(Float, Breaks[1])
	end

	def test_hyp_snps
		assert_kind_of(Array, Hyp)
		assert_kind_of(Float, Hyp[0])
	end

	def test_qq_cor
		hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
		fratio_breaks = SNPdist::fratio(hm, ht, Div, Genome_length)
		hyp = SNPdist::hyp_snps(fratio_breaks, Div, Genome_length)
		cor = SNPdist::qq_cor(Hyp, hyp)
		assert_kind_of(Float, cor)
		assert(0 <= cor && cor <= 1)
	end

end