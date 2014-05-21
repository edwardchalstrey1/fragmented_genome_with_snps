#encoding: utf-8
require_relative '../lib/write_it'
require 'test/unit'

class TestWriteIt < Test::Unit::TestCase

	File = "test/test/ratio_values.txt"

	def test_prime?
		assert_equal(true, WriteIt::prime?(7))
		assert_equal(false, WriteIt::prime?(4))
		assert_equal(true, WriteIt::prime?(2))
		assert_equal(true, WriteIt::prime?(971))
		assert_equal(false, WriteIt::prime?(972))
	end

	def test_file_to_array
		contents = WriteIt.file_to_array(File)
		assert_kind_of(Array, contents)
		assert_kind_of(String, contents[0])
	end

	def test_file_to_ints_array
		contents = WriteIt.file_to_ints_array(File)
		assert_kind_of(Array, contents)
		assert_kind_of(Integer, contents[0])
	end

	def test_file_to_floats_array
		contents = WriteIt.file_to_floats_array(File)
		assert_kind_of(Array, contents)
		assert_kind_of(Float, contents[0])
	end

end

