#encoding: utf-8

class FitnessScore
	require 'rinruby'

	# Input 0: Array of SNP positions
	# Input 1: Number of breaks (divisions) at which to bin intervals
	# Input 2: The length of the genome
	# Output 0: Array of number of SNPs in each bin
	# Output 1: Array of the breaks in the genome (the intervals created by div)
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

	# Input 0: Array where each value is the no. of homozygous SNPs in a division of the genome
	# Input 1: Array where each value is the no. of heterozygous SNPs in a division of the genome
	# Output: Array of ratios (floats) of homozygous to heterozygous SNPs for each division/bin of the genome
	def self.ratio(hm_count, ht_count)
		x = 0
		ratios = []
		hm_count.length.times do # the number of divisions of the genome (div)
			if hm_count[x] == 0 || ht_count[x] == 0
				ratios << 'NaN' # cannot calculate ratios
			else
				ratios << hm_count[x].to_f / ht_count[x].to_f
			end
			x+=1
		end
		return ratios
	end

	# Input 0: Array of expected ratios (floats) of homozygous to heterozygous SNPs for each division/bin of the genome
	# Input 1: Array of measured ratios (floats) of homozygous to heterozygous SNPs for each division/bin of the permutation
	# Output: Float between 0.0 and 1.0 where closely matching inputs are closer to 1.0 (pearson correlation)
	def self.score(expected, permutation)
		puts "expected: #{expected}"
		puts "permutation: #{expected}"
		x, ex, perm = 0, expected.dup, permutation.dup
		ex.length.times do
			if ex[x] == 'NaN' || perm[x] == 'NaN'
				ex.delete_at(x); perm.delete_at(x)
			else
				x+=1
			end
		end
		puts "ex: #{ex}"
		puts "perm: #{perm}"
		myr = RinRuby.new(echo = false)
		myr.assign 'x', ex
		myr.assign 'y', perm
		myr.eval 'score <- abs(cor(x,y))'
		fitness_score = myr.pull 'score'
		myr.quit
		# puts "fitness: #{fitness_score}"
		return fitness_score
	end
end