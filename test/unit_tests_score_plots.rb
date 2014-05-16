#encoding: utf-8
require_relative '../lib/score_plots/score_plots'
require 'test/unit'

class TestMetricPlot < Test::Unit::TestCase

	All_perms = MetricPlot::get_perms(2, 0, 1, 'small_dataset2', 'run13')

	def test_original_order
		ids = MetricPlot::original_order('small_dataset2')
		assert_kind_of(Array, ids)
		assert_kind_of(String, ids[0])
		assert_equal(53, ids.length)
	end

	def test_get_perms
		assert_equal(2, All_perms.length)
		assert_kind_of(Array, All_perms[0], 'population not array')
		assert_equal(54, All_perms[0][0].length, 'permutation wrong length') # frags plus fitness
		assert_kind_of(String, All_perms[0][0][0])
	end

	def test_plot_info
		x, y, se, group, best_sc = MetricPlot::plot_info(0, 1, 'lcs', All_perms, 'small_dataset2', 'run13')
		[x, y, se, group, best_sc].each do |i|
			assert_kind_of(Array, i)
			assert_equal(4, i.length)
		end
		[x, y, se, best_sc].each do |i|
			assert_in_delta(0.5, i[0], 0.5)
		end
	end
end