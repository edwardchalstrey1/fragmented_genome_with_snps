#encoding: utf-8
class RearrangementScore
	# Input 0: Array of objects in the correct order
	# Input 1: Array of the same objects as an incorrect permutation
	# Output 0: Ordinal similarity score value (0 = correct order)
	def self.rearrangement_score(original, permutation)
		both, difference_abs = [], []
		both << original.map{|x| permutation.index(x)} #works out the index of original values in permutation
		both << Array(0..(original.length - 1)) # index values that original originally at
		difference = both.transpose.map {|x| x.reduce(:-)} # taking away old position from new position, to find the distance that the frag has moved when re-ordered
		difference.each {|i| difference_abs << i.abs }
		return difference_abs.inject(:+) #high score = bad, score of 0 means the fragments in the right order
	end

	# Input 0: Array of objects in the correct order
	# Input 1: Array of the same objects as an incorrect permutation
	# Output 0: Ordinal similarity score value / worst possible score
	def self.score(original, permutation)
		"#{rearrangement_score(original, permutation)}/#{rearrangement_score(original, original.reverse)}"
	end

	def self.metric_2(original, permutation, version)
		pos_1 = permutation.index(original[0]) # position of the first object of original order in permutation
		if original[0] != permutation[0]
			new_a = [permutation[pos_1..-1], permutation[0..pos_1-1]].flatten # re-order the permutation to get first object at front
		else
			new_a = permutation
		end
		x = 0
		hds = [] # hamming distances
		new_a.each do |frag_id|
			if frag_id == original[x]
				hds << 0
			else
				hds << 1
			end
			x+=1
		end
		if version == 'b'
			score = (hds.inject(:+).to_f/permutation.length.to_f) + pos_1
		elsif version == 'c'
			score = hds.inject(:+) + (pos_1.to_f/permutation.length.to_f)
		else
			score = hds.inject(:+) + pos_1 # no need to take away 1 as the index in ruby is position - 1
		end
		return score
	end

	def self.metric_1_2_av(original, permutation)
		(rearrangement_score(original, permutation).to_f + metric_2(original, permutation, 'a').to_f)/2
	end

end