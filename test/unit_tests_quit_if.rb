#encoding: utf-8
require 'test/unit'
require_relative '../lib/quit_if'

class TestQuitIf < Test::Unit::TestCase

	def test_quit
		auc = QuitIf.quit([0.9,0.9,0.92,0.95,0.955])
		assert_kind_of(Float, auc)
		assert_equal(3.6975, ('%.4f' % auc).to_f)
	end
end