#encoding: utf-8

class LocateMutation
	require 'rinruby'

	# Input 0: List of SNP positions
	# Input 1: The number of equally spaced points at which the density is to be estimated. Specify n as a power of two.
	def self.find_peak(snps, n)
		myr = RinRuby.new(echo=false)
		myr.n = n
		myr.snps = snps
		myr.eval 'library(pracma)'
		myr.eval 'kernel_density <- density(snps, n=n)'
		myr.eval 'index <- match(max(kernel_density$y),kernel_density$y)'
		myr.eval 'peak <- kernel_density$x[index]'
		peak = myr.pull 'peak'
		myr.quit
		return peak.round # round to closest base position
	end
end