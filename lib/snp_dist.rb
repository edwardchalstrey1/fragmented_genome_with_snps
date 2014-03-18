#encoding: utf-8

class SNPdist
	require 'rinruby'
	
	### Model Genome ###

	# Make a list that models homozygous SNP positions
	def self.hm
		myr = RinRuby.new(echo = false)
		myr.eval 'hm <- rnorm(35, 10000000, 5000000)
		remove <- c()
		x <- 1
		for(i in hm){
		  if(i < 0 || i > 18585056){
		    remove <- c(remove, x)
		  }
		  x <- x+1
		}
		hm <- hm[-c(remove)]'
		hm = myr.pull 'hm'
		x = 0
		hm.each do |snp|
			hm[x] = snp.to_i
			x+=1
		end
		myr.quit
		return hm
	end

	# Make a list that models heterozygous SNP positions
	def self.ht
		myr = RinRuby.new(echo = false)
		myr.eval 'ht <- runif(3000, 1, 18585056)'
		ht = myr.pull 'ht'
		myr.quit
		return ht
	end

	### Fitness of re-order contigs algorithm ###

	# Gets the frequency of SNPs at breaks of div length for the homozygous and heterozygous SNP arrays
	# Divides to get the frequency ratio for each break (fratio)
	def self.fratio(hm, ht, div)
		myr = RinRuby.new(echo = false)
		myr.assign 'hm', hm
		myr.assign 'ht', ht
		myr.assign 'div', div
		myr.eval 'l <- 18585056
		breaks <- c(0)
		for(i in 1:div){
		  breaks <- c(breaks,(l/div)*i)
		}
		hmc <- hist(hm, breaks=breaks, plot=FALSE)$counts
		htc <- hist(ht, breaks=breaks, plot=FALSE)$counts
		'
		breaks = myr.pull 'breaks'
		hmc = myr.pull 'hmc'
		htc = myr.pull 'htc'
		myr.quit
		x = 0
		fratio = []
		hmc.length.times do
			count = (((hmc[x] + 1).to_f / (htc[x] + 1).to_f)*10).to_i # a measure of frequency
			if count == 0
				fratio << 0
			else
				fratio << count - 1
			end
			x+=1
		end
		return fratio, breaks
	end

	# Make hypothetical SNP array to test distribution of ratio
	def self.hyp_snps(fratio_breaks, div)
		hyp = []
		x = 0
		fratio_breaks[0].each do |freq|
			freq.times do
				hyp << rand(18585056.0/div.to_f) + fratio_breaks[1][x] # random value from within the range that the freq has been taken
			end
			x+=1
		end
		return hyp
	end

	# def self.qq(hyp)
	# 	myr = RinRuby.new(echo = false)
	# 	myr.assign 'hyp', hyp
	# 	myr.eval 'qqp <- qqnorm(hyp, plot.it=FALSE)'
	# 	myr.eval 'corr <- cor(qqp$x,qqp$y)'
	# 	corr = myr.pull 'corr'
	# 	myr.quit
	# 	return corr
	# end

	# Perform Shapiro-Wilks test for normality on a list
	def self.swtest(hyp)
		myr = RinRuby.new(echo = false)
		myr.assign 'hyp', hyp
		myr.eval 'test <- shapiro.test(hyp)'
		p_val = myr.pull 'test$p'
		myr.quit
		return p_val
	end

	# Make plots of the hypothetical SNP density and the ratio of hm to ht density to compare
	def self.plot_hyp(hyp, hm, ht)
		myr = RinRuby.new(echo = false)
		myr.assign 'hyp', hyp
		myr.assign 'hm', hm
		myr.assign 'ht', ht

		myr.eval 'png("~/fragmented_genome_with_snps/test/hypothetical_snps/hyp.png")
		plot((1:512)*36298.9375, density(hyp)$y, xlab="Hyp snps dist")
		dev.off()'

		myr.eval 'hmd <- density(hm, from=0, to=18585056)
		htd <- density(ht, from=0, to=18585056)
		ratio <- hmd$y/htd$y'

		myr.eval 'png("~/fragmented_genome_with_snps/test/hypothetical_snps/d_ratio.png")
		plot((1:512)*36298.9375, ratio)
		dev.off()'
		myr.quit
	end
end

hm = SNPdist::hm
ht = SNPdist::ht

fratio_breaks = SNPdist::fratio(hm, ht, 1000)
# puts fratio_breaks[0]
hyp = SNPdist::hyp_snps(fratio_breaks, 1000)
puts hyp == hyp.uniq
puts hyp.length
# puts hyp

puts SNPdist::swtest(ht)
puts SNPdist::swtest(hm)
puts SNPdist::swtest(hyp)

SNPdist::plot_hyp(hyp,hm,ht)