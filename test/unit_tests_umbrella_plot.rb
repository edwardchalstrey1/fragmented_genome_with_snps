#encoding: utf-8
require_relative '../lib/score_plots/umbrella_plot.rb'
require 'test/unit'

class TestUPlot < Test::Unit::TestCase

	def setup
		@dataset = 'small_dataset2a'
		@runs = UPlot.get_runs(@dataset)
	end
	
	def test_get_runs
		assert_equal(4, @runs.length)
		assert_equal('p_run1', @runs[0])
	end

	def test_get_gens
		gens = UPlot.get_gens(@dataset, @runs[0])
		assert_equal(6, gens)
	end

	def test_plot_info
		gens, fits, all_runs, perms, param_types = UPlot.plot_info(@dataset)

		assert_equal(11060, gens.length) # Total permutations in plot
		assert_equal(11060, fits.length)
		assert_equal(11060, all_runs.length)
		assert_equal(11060, perms.length)
		assert_equal(11060, param_types.length)

		assert_equal(0, gens[0])
		assert_equal(0.9521, ('%.4f' % fits[0]).to_f) # First permutation fitness
		assert_equal('p_run1', all_runs[0])
		assert_equal(53, perms[0].length)
		assert_equal('frag35', perms[0][0])
		assert_equal('p1', param_types[0])

		assert_equal(199, gens[-1])
		assert_equal(0.6832, ('%.4f' % fits[-1]).to_f) # Permutation 9 fitness, p_run 50: permutation files ordered alphabetically (order doesn't matter)
		assert_equal('p_run50', all_runs[-1])
		assert_equal(53, perms[-1].length)
		assert_equal('frag17', perms[-1][0])
		assert_equal('p5', param_types[-1])
	end

end