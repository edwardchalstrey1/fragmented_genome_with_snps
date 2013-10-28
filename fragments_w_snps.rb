require "rubygems"
require "rinruby"
require "json"

def generate_positions (breaks, counts)
	snp_pos = [] #each of the ranges between each of our 'breaks' (between bars of the histogram) needs a proportion of snps
	x = Random.new
	i = 0
	counts.each do |l| #for each count do
		if l != nil
			(breaks[i].to_i..breaks[i+1].to_i).each do |j| #for each range we need a random no. for each position (ranges in this case are 10k)
				newr = x.rand(1..200000) #the upper bound is the length of the sequence, where l is the frequency (if we multiply each frequency by the same constant...)
				if newr <= l*50
					snp_pos << j #add the position (j) to array if its random no. (newr) is <= count. 
					#count is frequency, so there should be count/200k snp positions in each range
				end
			end
		end
		i+=1
	end
	return snp_pos
end
def sample (snp_pos, rinruby)
	rinruby.assign "pos_l", snp_pos #take 1000 snps (still normally distributed) to use
	rinruby.eval "pos_s <- sample(pos_l, 1000)"
	final_positions = rinruby.pull "pos_s"
	return final_positions
end
def normal_dist
	myr = RinRuby.new(echo = false)
	myr.eval "x <- rnorm(200000, 100000, 19000)"
	myr.eval "h <- hist(x, breaks = 8000, plot = FALSE)" # The larger the number of breaks, the closer to a normal distribution the snp positions will be see generate positions^
	myr.eval "b <- h$breaks"
	myr.eval "c <- h$counts"
	c = myr.pull "c"
	b = myr.pull "b"
	snp_pos = generate_positions(b, c)
	puts "The 1000 SNPs are sampled from "+snp_pos.length.to_s+" normally distributed positions"
	s = sample(snp_pos, myr)
	return s
end

def unique? (snp_pos) #This method is an experiment/test to show that the previous methods have produced unique snp positions
	unique_pos = []
	snp_pos.each do |q|
		if unique_pos.include? q
			next
		else
			unique_pos << q
		end
	end
	puts "All snp positions generated are unique: " + (snp_pos.length == unique_pos.length).to_s
end

def make_snp_seq (snp_pos)
	snp_seq = ['a', 'g', 't', 'g', 'c']*40000 #This step is actually unecessary, though it shows later on when creating a vcf file that all works.
	return snp_seq
end
def get_frags (seq)
	frags = []
	all_frags_pos = []
	rt = 0
	while rt < seq.length
		frag_length = rand(200) + 50
		frag = seq[rt..(rt+frag_length)]
		rt = rt+frag_length
		frags << frag
	end
	return frags
end
def remove_extra_frags (frags)
	x = frags.flatten.length - 200000
	y = frags[-1].length
	c = 2
	while y < x
		y = y + frags[-c].length
		c+=1
	end
	frags2 = frags - frags[-c..-1]
	leftover = frags.flatten.length - frags2.flatten.length
	leftover2 = frags.flatten.length - 200000
	final_frag = ['a']*(leftover - leftover2)
	frags2 << final_frag
	return frags2
end
def pos_each_frag (snp_pos, frags) #get the positions for snps on individual frags
	sorted_pos = snp_pos.sort
	p_ranges = []
	frags.each do |i|
		p_ranges << i.length + p_ranges[-1].to_i #adding the fragment lengths to get the upper bounds of ranges of positions on the original seq.
	end
	first_pos = [] #then, to work out the first position of each fragment
	first_pos << 0
	first_pos << p_ranges[0..-2]
	first_pos = first_pos.flatten
	p = 0
	t = 0
	all_frags_pos = [] # the positions of snps on each fragment, array of arrays
	p_ranges.each do |jj| #for each of the upper ranges (j) that the positions could be in
		each_frag_pos = []
		while sorted_pos[t].to_i < jj && !sorted_pos[t].nil? do #make the loop quit before trying to access index of snp_pos that doesn't exist
			each_frag_pos << sorted_pos[t].to_i 	# add all of the positions < jj to a new array for each frag
			t += 1
		end
		z = 0
		y = first_pos[p].to_i
		each_frag_pos.each do |e|   
			each_frag_pos[z] = e - y #taking away the value at the first position of each fragment to get the true positions
			z += 1
		end
		p += 1
		all_frags_pos << each_frag_pos # add the positions for the frag to an array of the each_frag arrays
	end
	return all_frags_pos
end
def storage (frags, pos)
	frags_with_positions = []
	x = 0
	frags.each do |i|
		frag_pos_hash = {}
		frag_pos_hash[:frag] = i.join.to_s
		frag_pos_hash[:pos] = pos[x]
		frags_with_positions << frag_pos_hash  #so ideally i would use the hashed pos_each_frag method above, but for some reason dnt work
		x+=1
	end
	return frags_with_positions
end
def write_json (array)
	File.open("frags_with_positions.json", "w") do |f|
		f.write(array.to_json)
	end
end

snp_pos = normal_dist #make a normal distribution of snp positions
unique?(snp_pos) #are the positions generated unique?

snp_seq = make_snp_seq(snp_pos) #make a sequence with snps at the positions^
frags_plus = get_frags(snp_seq) #get an array of fragments

frags = remove_extra_frags(frags_plus)

positions_each = pos_each_frag(snp_pos, frags)

frags_with_positions = storage(frags, positions_each)
write_json(frags_with_positions)

File.open("snp_pos.txt", "w+") do |f|
	snp_pos.each { |i| f.puts(i) }
end

