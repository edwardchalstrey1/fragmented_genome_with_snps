#encoding: utf-8

class SNPdist
	require 'rinruby'

	# Gets the frequency of SNPs at breaks of div length for the homozygous and heterozygous SNP arrays
	# Divides to get the frequency ratio for each break (fratio), div = no. of breaks
	def self.fratio(hm, ht, div, genome_length)
		myr = RinRuby.new(echo = false)
		myr.assign 'hm', hm
		myr.assign 'ht', ht
		myr.assign 'div', div
		myr.assign 'l', genome_length
		myr.eval 'breaks <- c(0)
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
	# div = number of breaks
	def self.hyp_snps(fratio_breaks, div, genome_length)
		hyp = []
		x = 0
		fratio_breaks[0].each do |freq| 
			freq.times do
				hyp << rand(genome_length/div.to_f) + fratio_breaks[1][x] # random value from within the range that the freq has been taken
			end
			x+=1
		end
		return hyp
	end

	# Make plots of the hypothetical SNP density and the ratio of hm to ht density to compare
	def self.plot_hyp(hyp, location, dataset, gen, genome_length)
		myr = RinRuby.new(echo = false)
		myr.assign 'hyp', hyp
		myr.assign 'location', location
		myr.assign 'dataset', dataset
		myr.assign 'gen', gen
		myr.assign 'genome_length', genome_length
		myr.eval 'png(paste("~/",location,"/", dataset,"/gen_", gen, "_best_permutation_distribution.png", sep=""))
		plot((1:512)*(genome_length/512), density(hyp)$y, xlab="Approximated ratio of homozygous to heterozygous SNP density")
		dev.off()'		
		myr.quit
	end

	# The old way of calculating density, should be similar looking plot to plot_hyp if hyp snps are good
	def self.plot_dens(hm, ht, location, dataset, genome_length)
		myr = RinRuby.new(echo = false)
		myr.assign 'hm', hm
		myr.assign 'ht', ht
		myr.assign 'location', location
		myr.assign 'dataset', dataset
		myr.assign 'genome_length', genome_length
		myr.eval 'hmd <- density(hm, from=0, to=genome_length)
		htd <- density(ht, from=0, to=genome_length)
		ratio <- hmd$y/htd$y'
		myr.eval 'png(paste("~/",location,"/", dataset,"/density_vector_ratio.png", sep=""))
		plot((1:512)*(genome_length/512), ratio, xlab="Ratio of homozygous to heterozygous SNP density")
		dev.off()'
		myr.quit
	end

	# Q-Q correlaion: ### TODO get rid of plot
	def self.qq_cor(example_ratio, ratio)
		myr = RinRuby.new(echo = false)
		myr.assign 'ex', example_ratio
		myr.assign 'ra', ratio
		myr.eval '# png("~/fragmented_genome_with_snps/test/hypothetical_snps/qqplot.png")
		qqp <- qqplot(ex, ra, plot.it=FALSE)
		# dev.off()
		corr <- cor(qqp$x,qqp$y)'
		corr = myr.pull 'corr'
		myr.quit
		return corr
	end
end