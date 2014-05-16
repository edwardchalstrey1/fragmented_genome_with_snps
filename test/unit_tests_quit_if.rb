#encoding: utf-8
require 'test/unit'
require_relative '../lib/quit_if'

class TestQuitIf < Test::Unit::TestCase

	def test_quit
		slope = QuitIf.quit([0.9,0.9,0.92,0.95,0.955])
		assert_kind_of(Float, slope)
		assert_in_delta(0.5, slope, 0.5)
	end
end