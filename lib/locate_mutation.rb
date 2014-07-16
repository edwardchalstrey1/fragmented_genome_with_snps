#encoding: utf-8

class LocateMutation
	require 'rinruby'

	def self.find_peak(snps)
		myr = RinRuby.new(echo=false)
		myr.snps = snps
		myr.eval 'library(pracma)'
		myr.eval 'kernel_density <- density(snps)'
		myr.eval 'index <- match(max(kernel_density$y),kernel_density$y)'
		myr.eval 'peak <- kernel_density$x[index]'
		peak = myr.pull 'peak'
		myr.quit
		return peak
	end
end