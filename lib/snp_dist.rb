#encoding: utf-8

class SNPdist
	require 'rinruby'

	### Hypothetical SNP positions ###

	# Input 0: Array of ratios for homozygous to heterozygous SNPs at divisions of the genome
	# Input 1: The length of the genome
	# Output: A list of "hypothetical" SNP positions which represents the distribution of homozygous/heterozygous SNP density ratio
	def self.hyp_snps(ratios, genome_length)
		breaks = []
		(1..ratios.length).to_a.each do |i|
			breaks << (genome_length/ratios.length.to_f)*i
		end
		hyp, x = [], 0
		ratios.each do |ratio| 
			(ratio*10).to_i.times do
				hyp << rand(genome_length/ratios.length.to_f) + breaks[x] # random value from within the range that the freq has been taken
			end
			x+=1
		end
		return hyp # These don't need to be unique or integers like the real SNPs, since they are just representing a distribution
	end

	### Plotting Methods ###

	# Input 0: Array of ratios for homozygous to heterozygous SNPs at divisions of the genome
	# Input 1: Location at which to save the plot
	# Input 2: The dataset to use (and the sub-folder to save the plot in)
	# Input 3: The generation of the genetic algorithm the ratio is being plotted for
	# Input 4: The length of the genome
	# Output: Plot: ratio of homozygous to heterozygous SNP density
	def self.plot_ratio(ratios, location, dataset_run, gen, genome_length)
		myr = RinRuby.new(echo = false)
		myr.assign 'ratios', ratios
		myr.assign 'location', location
		myr.assign 'dataset_run', dataset_run
		myr.assign 'gen', gen
		myr.assign 'genome_length', genome_length
		myr.eval 'png(paste("~/",location,"/", dataset_run,"/Gen", gen, "_lists/best_permutation_ratios_", (length(ratios)/1000), "Kdiv.png", sep=""))
		plot((1:length(ratios))*(genome_length/length(ratios)), ratios, xlab=paste("Genome (contigs ordered by best permutation in generation ", gen, ")", sep=""), 
			ylab="Ratio",
			main=paste("Ratio of homozygous to heterozygous SNP density
			 calculated at ", length(ratios), " divisions of the genome", sep=""))
		dev.off()'
		myr.quit
	end

	# Input 0: A list of SNP positions 
	# Input 1: Location at which to save the plot
	# Input 2: The dataset to use (and the sub-folder to save the plot in)
	# Input 3: The generation of the genetic algorithm the ratio is being plotted for
	# Input 4: The length of the genome
	# Input 5: String indicating the type of SNPs that Input 0 are
	# Input 6: Title of plot 
	# Output: Plot of kernel density estimate for the SNPs over the genome
	def self.plot_snps(snp_pos, location, dataset_run, gen, genome_length, type, title)
		myr = RinRuby.new(echo = false)
		myr.assign 'snp_pos', snp_pos
		myr.assign 'location', location
		myr.assign 'dataset_run', dataset_run
		myr.assign 'gen', gen
		myr.assign 'genome_length', genome_length
		myr.assign 'type', type
		myr.assign 'title', title
		myr.eval 'png(paste("~/",location,"/", dataset_run,"/Gen", gen, "_lists/best_permutation_distribution_", type, ".png", sep=""))
		plot((1:512)*(genome_length/512), density(snp_pos)$y, xlab=paste("Genome (contigs ordered by best permutation in generation ", gen, ")", sep=""),
			ylab="Density", main=title)
		dev.off()'		
		myr.quit
	end
	
	# Input 0: List of homozygous SNP positions
	# Input 1: List of heterozygous SNP positions
	# Input 2: Location at which to save the plot
	# Input 3: The dataset to use (and the sub-folder to save the plot in)
	# Input 4: The length of the genome
	# Input 5: Name of plot
	# Ootput: Plot of ratio of hom/het SNP density
	def self.plot_dens(hm, ht, location, dataset, genome_length, name)
		myr = RinRuby.new(echo = false)
		myr.assign 'hm', hm
		myr.assign 'ht', ht
		myr.assign 'location', location
		myr.assign 'dataset', dataset
		myr.assign 'genome_length', genome_length
		myr.assign 'name', name
		myr.eval 'hmd <- density(hm, from=0, to=genome_length)
		htd <- density(ht, from=0, to=genome_length)
		ratio <- hmd$y/htd$y'
		myr.eval 'png(paste("~/",location,"/", dataset,"/", name, ".png", sep=""))
		plot((1:512)*(genome_length/512), ratio, xlab="Genome", ylab="Ratio of homozygous to heterozygous SNP density")
		dev.off()'
		myr.quit
	end
end