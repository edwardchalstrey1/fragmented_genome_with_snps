original_order = [1,2,3,4,5,6,7,8,9,10]
#rearranged = frags_original_order.shuffle
#rearranged = frags_original_order.reverse
rearranged = [9,10,8,7,6,5,4,3,1,2]

















position = original_order.map{|x| rearranged.index(x)} #works out the index of original_order values in rearranged
index_values = Array(0..(original_order.length - 1)) # index values that fasta_ids originally at
both = []
both << position
both << index_values
difference = both.transpose.map {|x| x.reduce(:-)} # taking away old position from new position, to find the distance that the frag has moved when re-ordered
difference_abs = []
difference.each do |i|
	difference_abs << i.abs
end
score = difference_abs.inject(:+)
puts score