#encoding: utf-8

class SNPdist
	require 'rinruby'
	
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

	def self.ht
		myr = RinRuby.new(echo = false)
		myr.eval 'ht <- runif(3000, 1, 18585056)'
		ht = myr.pull 'ht'
		myr.quit
		return ht
	end

	def self.fratio(hm, ht)
		myr = RinRuby.new(echo = false)
		myr.assign 'hm', hm
		myr.assign 'ht', ht
		myr.eval 'l <- 18585056
		div <- 100 #l/10000
		breaks <- c(0)
		for(i in 1:div){
		  breaks <- c(breaks,(l/div)*i)
		}
		fhm <- summary(cut(hm, breaks))
		fht <- summary(cut(ht, breaks))
		fratio <- fhm/fht'
		fratio = myr.pull 'fratio'
		return fratio
	end

	# Make hypothetical SNP array to test distribution of ratio
	def self.hyp_snps(fratio)
		hyp = []
		fratio.each do |freq|
			freq.times do
				hyp << # random value from within the range that the freq has been taken
			end
		end
		return hyp
	end

end

hm = SNPdist::hm
ht = SNPdist::ht

puts SNPdist::fratio(hm,ht)