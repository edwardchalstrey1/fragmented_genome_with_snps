#encoding: utf-8
class RearrangementScore
	# Input 0: Array of objects in the correct order
	# Input 1: Array of the same objects as an incorrect permutation
	# Output 0: Ordinal similarity score value (0 = correct order)
	def self.rearrangement_score(frags_original_order, rearranged)
		position_each_frag_id_in_d = frags_original_order.map{|x| rearranged.index(x)} #works out the index of frags_original_order values in rearranged
		index_values = Array(0..(frags_original_order.length - 1)) # index values that frags_original_order originally at
		both = []
		both << position_each_frag_id_in_d
		both << index_values
		difference = both.transpose.map {|x| x.reduce(:-)} # taking away old position from new position, to find the distance that the frag has moved when re-ordered
		difference_abs = []
		difference.each do |i|
			difference_abs << i.abs
		end
		score = difference_abs.inject(:+) #high score = bad, score of 0 means the fragments in the right order
		return score
	end

	def self.score(frags_original_order, rearranged)
		score = rearrangement_score(frags_original_order, rearranged)
		worst = rearrangement_score(frags_original_order, frags_original_order.reverse)
		return "#{score}/#{worst}"
	end
end