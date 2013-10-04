require "rubygems"
require "rinruby"
require "json"

def generate_positions (breaks, counts)
	snp_pos = [] #each of the ranges between each of our 'breaks' (between bars of the histogram) needs a proportion of snps
	x = Random.new
	i = 0
	counts.each do |l| #for each count do
		(breaks[i].to_i..breaks[i+1].to_i).each do |j| #for each range we need a random no. for each position (ranges in this case are 10k)
			newr = x.rand(1..200000) #make a random no. within range 
			if newr <= l
				snp_pos << j #add the position (j) to array if its random no. (newr) is <= count. 
				#count is frequency, so there should be count/200k snp positions in each range
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
	myr.eval "h <- hist(x, plot = FALSE)"
	myr.eval "b <- h$breaks"
	myr.eval "c <- h$counts"
	c = myr.pull "c"
	b = myr.pull "b"
	#puts snp_counts.length # 18 counts because split into 20 lots of 10kb for 200kb sequence, 2 counts have 0 snps (one at either end of dist.)
	#puts snp_breaks.length # 19 breaks because these are the breaks between the 20 'bars' of the histogram
	snp_pos = generate_positions(b, c)
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
	a = ['a', 'g', 't', 'g', 'a'] #This step is actually unecessary, though it shows later on when creating a vcf file that all works.
	snp_seq = a*200000
	snp_pos.each do |i|
		snp_seq[i] = 'c' #replacing the slots at the snp positions with snp's (always c for simplicity)
	end
	return snp_seq
end
def get_frags (seq, snp_pos)
	frags = []
	all_frags_pos = []
	rt = 0
	while rt < seq.length
		frag_length = rand(200) + 50
		frag = seq[rt..(rt+frag_length)]
		pos = snp_pos.select{|x| x>rt and x<(rt+frag_length)}
		rt = rt+frag_length
		frags << frag
		all_frags_pos << pos
	end
	return frags
end

snp_pos = normal_dist
unique_test = unique?(snp_pos)
y = make_snp_seq(snp_pos)
x = get_frags(y, snp_pos)
puts x.flatten.length






