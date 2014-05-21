#encoding: utf-8
require 'rinruby'

class QuitIf

	# Input: Array of number values
	# Output: The slope coefficient of the input array plotted against 1..input_length
	def self.quit(best_permutations)
		myr = RinRuby.new(echo = false)
		myr.assign 'y', best_permutations
		myr.eval 'x <- 1:length(y)'
		myr.eval 'slope <- lm(y~x)$coefficients[[2]]'
		slope = myr.pull 'slope'
		return slope
	end
end