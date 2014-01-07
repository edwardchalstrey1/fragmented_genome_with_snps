def division (frags, x, prime) #number of frags
	if prime == "n"
		until frags.length/x == (frags.length/x).to_i && x == x.to_i && x <= frags.length
			x = (frags.length/10).to_f + rand(frags.length).to_f
		end
	elsif prime == "p"
		until frags.length/x == (frags.length/x).to_i && x == x.to_i && x < frags.length
			x = (frags.length/10).to_f + rand(frags.length).to_f
		end
	end		
	return x
end

m = 10#1303
mum = (1..m).to_a
dad = (m+1..m*2).to_a
x = division(mum, 1.5, "n")

if x == mum.length # to compensate for datasets with a prime number of fragments:
	ig = rand(mum.length)-1 # choose a random element of the fasta array to ignore
	mig = mum[ig] # we can add these frags back at their original positions after recombination
	dig = dad[ig]
	mum.delete_at(ig)
	dad.delete_at(ig)
	x = division(mum, 1.5, "p")
end

puts x
puts mum.length/x

mum = mum.each_slice(x).to_a
dad = dad.each_slice(x).to_a
puts mum.length
puts dad[0].length
puts ig != nil

# Let's say we use one chunk of the dad solution, the rest mum

ch = rand(dad.length)-1 # choose one of the chunks of fragments to keep from dad
child = mum.dup.flatten
y = 0
mum[ch].each do |frag| # place each frag in the equivalent mum chunk into the position it's corresponding frag (from dad) occupies in mum
	chunk_frag = dad[ch][y]
	pos = mum.flatten.index(chunk_frag) # MAKE THIS WORK!!!
	child[pos] = frag
end