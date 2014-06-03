#encoding: utf-8

class QuitIf
	require 'rinruby'
	# Input: Array of number values
	# Output: The area under the curve for the input array plotted against 1..input_length
	def self.quit(fitness_scores)
		myr = RinRuby.new(echo = false)
		myr.assign 'y', fitness_scores
		myr.eval 'x <- 1:length(y)'
		myr.eval 'library("pracma")'
		myr.eval 'require(pracma)'
		myr.eval 'auc <- trapz(x,y)'
		auc = myr.pull 'auc'
		myr.quit
		return auc
	end
end