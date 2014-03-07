#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/write_it'
require 'test/unit'

class TestWriteIt < Test::Unit::TestCase

	def test_prime?
		assert_equal(true, WriteIt::prime?(7))
		assert_equal(false, WriteIt::prime?(4))
		assert_equal(true, WriteIt::prime?(2))
		assert_equal(true, WriteIt::prime?(971))
		assert_equal(false, WriteIt::prime?(972))
	end

end

