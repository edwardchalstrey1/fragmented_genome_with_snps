require "rubygems"
require "rinruby"

def skew_grad (min_snps_per_frag)
	myr = RinRuby.new(echo = false)
	myr.assign "min", min_snps_per_frag
	myr.eval "source('fragmented_genome_with_snps/skew_scatta.R')"
	myr.eval "gradients <- skew_grad(min)"
	myr.eval "gradients <- as.numeric(gradients)[1]"
	myr.eval "x <- 2.043065e-09"
	return myr.pull "gradients"
end

puts skew_grad(6)