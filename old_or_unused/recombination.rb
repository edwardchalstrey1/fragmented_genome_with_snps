def prime? (n)
	for d in 2..(n - 1)
		if (n % d) == 0
    		return false
    	end
	end
	true
end

def division (frags, prime) #number of frags
	x = 1.5
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

def recombine (mum, dad)
	kid = []
	1.times do
		mum1 = "bogey"
		x = division(mum, "n")
		if x == mum.length && prime?(x) == false
			puts "redoing1..."
			redo
		elsif x == mum.length # to compensate for datasets with a prime number of fragments:
			puts "we gotta prime"
			ig = rand(mum.length)-1 # choose a random element of the fasta array to ignore
			mig = mum[ig] # we can add these frags back at their original positions after recombination
			dig = dad[ig]
			mum1 = mum.dup
			dad1 = dad.dup
			mum1.delete_at(ig)
			dad1.delete_at(ig)
			x = division(mum1, "p")
		end
		if mum1 != "bogey"
			mum1 = mum1.each_slice(x).to_a
			dad1 = dad1.each_slice(x).to_a
		else
			mum1 = mum.each_slice(x).to_a
			dad1 = dad.each_slice(x).to_a
		end
		puts "Ignored a frag? " +(ig != nil).to_s
		# Let's say we use one chunk of the dad solution, the rest mum
		puts "x is " +x.to_s
		ch = rand(dad1.length)-1 # choose one of the chunks of fragments to keep from dad
		child = mum1.dup.flatten
		y = 0
		pos_array = []
		mum1[ch].each do |frag| # place each frag in the equivalent mum chunk into the position it's corresponding frag (from dad) occupies in mum
			chunk_frag = dad1[ch][y] # the equivalent frag in the chosen dad chunk
			pos = mum1.flatten.index(chunk_frag) # the position of the dad chunk frag in mum
			c_pos = mum1.flatten.index(frag) # the position of the frag in mum
			pos_array << pos
			#puts pos
			y+=1
		end
		if pos_array.include?(nil)
			puts "redoing2..."
			redo
		else
			y = 0
			pos_array.each do |pos|
				unless dad1[ch].include?(mum1[ch][y])
					child[pos] = mum1[ch][y]
					child[mum1.flatten.index(mum1[ch][y])] = dad1[ch][y] # swapping the positions of the frag and chunk frag, to give their positions in child
				end
				y+=1
			end
		end
		puts "kid length same: "+(mum1.flatten.length == child.length).to_s
		if ig != nil
			if dad1[ch].include?(dig) # add the ignored fragment at it's 
				child.insert(ig, dig)
			else
				child.insert(ig, mig)
			end
		end
		puts "kid unique? "+(child == child.uniq).to_s
		kid << child
	end
	return kid[0]
end

m = 1303
mum = (1..m).to_a
dad = mum.dup.shuffle
recombine(mum, dad)