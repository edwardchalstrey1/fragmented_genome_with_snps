#encoding: utf-8
require_relative '../lib/score_plots/umbrella_plot.rb'
require 'test/unit'

class TestUPlot < Test::Unit::TestCase

	def setup
		@dataset = 'small_dataset2a'
		@runs = UPlot.get_runs(@dataset)
	end
	
	def test_get_runs
		assert_equal(2, @runs.length)
		assert_equal('p_run1', @runs[0])
	end

	def test_get_gens
		gens = UPlot.get_gens(@dataset, @runs[0])
		assert_equal(6, gens)
	end

	def test_plot_info
		gens, fits, runs = UPlot.plot_info(@dataset)

		assert_equal(560, gens.length) # Total permutations in plot
		assert_equal(560, fits.length)
		assert_equal(560, runs.length)

		assert_equal(0, gens[0])
		assert_equal(0.0479, ('%.4f' % fits[0]).to_f) # First permutation fitness
		assert_equal('p_run1', runs[0])

		assert_equal(24, gens[-1])
		assert_equal(0.6264, ('%.4f' % fits[-1]).to_f) # Permutation 9 fitness, permutation files ordered alphabetically (order doesn't matter)
		assert_equal('p_run2', runs[-1])
	end

end