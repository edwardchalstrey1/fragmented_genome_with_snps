#encoding: utf-8

class FitnessScore
	require 'rinruby'

	# Input 0: Array of SNP positions
	# Input 1: Number of breaks (divisions) in the genome to count the number of SNPs in
	# Input 2: The length of the genome
	# Output: Array of number of SNPs in each genome division
	def self.count(snp_pos, div, genome_length)
		myr = RinRuby.new(echo = false)
		myr.assign 'snp_pos', snp_pos
		myr.assign 'div', div
		myr.assign 'l', genome_length
		myr.eval 'breaks <- c(0)
		for(i in 1:div){
		  breaks <- c(breaks,(l/div)*i)
		}
		counts <- hist(snp_pos, breaks=breaks, plot=FALSE)$counts'
		counts = myr.pull 'counts'
		myr.quit
		return counts
	end

	# Input 0: Array of homozygous SNP positions
	# Input 1: Array of heterozygous SNP positions
	# Input 2: Number of breaks (divisions) in the genome to count the number of SNPs in
	# Input 3: The length of the genome
	# Output: Array of ratios (floats) of homozygous to heterozygous SNPs for each division of the genome,
	#         where 1 has been added to each count, to avoid ratios of infinity or zero
	def self.ratio(hm, ht, div, genome_length)
		hm_count = FitnessScore::count(hm, div, genome_length)
		ht_count = FitnessScore::count(ht, div, genome_length)
		x = 0
		ratios = []
		hm_count.length.times do # the number of divisions of the genome (div)
			count_ratio = ((hm_count[x] + 1).to_f / (ht_count[x] + 1).to_f) # a measure of ratio
			ratios << count_ratio
			x+=1
		end
		return ratios
	end

	# Input 0: Array of expected ratios (floats) of homozygous to heterozygous SNPs for each division/bin of the genome
	# Input 1: Array of measured ratios (floats) of homozygous to heterozygous SNPs for each division/bin of the permutation
	# Output: Float between 0.0 and 1.0 where closely matching inputs are closer to 1.0 (pearson correlation)
	def self.score(expected, permutation)
		myr = RinRuby.new(echo = false)
		myr.assign 'x', expected
		myr.assign 'y', permutation
		myr.eval 'score <- abs(cor(x,y))'
		fitness_score = myr.pull 'score'
		myr.quit
		return fitness_score
	end

	def self.distance_score(hm)
		score = 0
		hm.each_cons(2).map { |a,b| score+=(b-a) }
		score.to_f
	end

end