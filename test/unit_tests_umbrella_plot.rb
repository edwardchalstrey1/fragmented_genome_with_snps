#encoding: utf-8
require_relative '../lib/score_plots/umbrella_plot.rb'
require 'test/unit'

class TestUPlot < Test::Unit::TestCase
	def setup
	end
	
	def test_get_dirs
		dirs = UPlot.get_dirs('10K_dataset3')
		assert_equal(4, dirs.length)
	end
end