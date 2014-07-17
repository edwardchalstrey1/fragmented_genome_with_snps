#encoding: utf-8
require_relative '../lib/locate_mutation'
require 'test/unit'

class TestLocateMutation < Test::Unit::TestCase

	def test_find_peak
		snps = [105,109,87,96,110,95,100,88,110,92]
		assert_equal(93.84763, ('%.5f' % LocateMutation.find_peak(snps, 512)).to_f)
	end

end