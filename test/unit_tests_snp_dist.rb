#encoding: utf-8
require_relative '../lib/snp_dist'
require_relative '../lib/write_it'
require 'test/unit'

class TestSNPdist < Test::Unit::TestCase

	def setup
		@ratios = WriteIt.file_to_floats_array("test/test/ratios_example.txt")
		@hyp = SNPdist.hyp_snps(@ratios, 2000)
	end

	def test_hyp_snps
		assert_kind_of(Array, @hyp)
		assert_kind_of(Float, @hyp[0])
	end

	def test_ylim
		ylim_ratio = SNPdist.get_ylim(@hyp, 2000, 'ratio')
		ylim_density = SNPdist.get_ylim(@hyp, 2000, 'density')
		assert_kind_of(Float, ylim_ratio)
		assert_kind_of(Float, ylim_density)
	end
end