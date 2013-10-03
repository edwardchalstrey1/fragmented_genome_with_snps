require "rubygems"
require "rinruby"
require "json"

myr = RinRuby.new(echo = false)
myr.eval "x <- rnorm(200000, 100000, 19000)"
myr.eval "h <- hist(x, plot = FALSE)"
myr.eval "b <- h$breaks"
myr.eval "c <- h$counts"
snp_counts = myr.pull "c"
snp_breaks = myr.pull "b"
#puts snp_counts.length # 18 counts because split into 20 lots of 10kb for 200kb sequence, 2 counts have 0 snps (one at either end of dist.)
#puts snp_breaks.length # 19 breaks because these are the breaks between the 20 'bars' of the histogram
snp_pos = [] #each of the ranges between each of our 'breaks' (between bars of the histogram) needs a proportion of snps
x = Random.new
i = 0
snp_counts.each do |l| #for each count do
	(snp_breaks[i].to_i..snp_breaks[i+1].to_i).each do |j| #for each range we need a random no. for each position (ranges in this case are 10k)
		newr = x.rand(1..200000) #make a random no. within range 
		if newr <= l
			snp_pos << j #add the position (j) to array if its random no. (newr) is <= count. 
			#count is frequency, so there should be count/200k snp positions in each range
		end
	end
	i+=1
end
myr.assign "pos_l", snp_pos #take 1000 snps (still normally distributed) to use
myr.eval "pos_s <- sample(pos_l, 1000)"
final_positions = myr.pull "pos_s"
unique_pos = []
final_positions.each do |q|
	if unique_pos.include? q
		next
	else
		unique_pos << q
	end
end
puts "There are " + final_positions.length.to_s + " SNP\'s"
puts unique_pos.length.to_s + " of the SNP positions are unique" #this proves each snp position is unique
puts
#Create sequence with SNPs at designated positions
a = ['a', 'g', 't', 'g', 'a'] #This step is actually unecessary, though it shows later on when creating a vcf file that all works.
snp_seq = a*200000
unique_pos.each{
	|i|
	snp_seq[i] = 'c' #replacing the slots at the snp positions with snp's (always c for simplicity)
}

#THE LOOP BELOW WORKS, BUT IS VERY VERY SLOW
#Split sequences into random sized fragments
frags = []
until frags.flatten.join.to_s.split(//).length >= 200000 do #once the frags array has all 200Kb from snp_seq
	frags << (snp_seq.each_slice(rand(200) + 50).to_a)[0] #give one frag of 50-200b to the frags array
	seq_string = snp_seq.join.to_s  #create a string from the sequence for slicing
	frags_string = frags[-1].join.to_s #create a string of the most recently added fragment to slice off
	seq_string.slice! frags_string #slice the fragment from the start of the string and... 
	snp_seq = seq_string.split(//) #modify the sequence so it no longer contains the removed frags, then the loop can repeat
	puts "Progress: " + ((frags.flatten.join.to_s.split(//).length)/2000).to_s + "%"
end
all = frags.flatten
x = all.length - 200000
(frags[-1])[-(x)..-1] = 'z'
frags[-1].delete_if {|i| i > 'y'} #shorten the final fragment that has extra nucleotides due to the until >= 200Kb loop
puts "Total nucleotides: " + (frags.flatten.join.to_s.split(//).length).to_s
puts "There are a total of " + frags.length.to_s + " fragments"

#Find the positions of snps on the fragments, based on the positions in the un-fragmented sequence
sorted_pos = unique_pos.sort
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
puts "The positions are separated into " + (all_frags_pos.length).to_s + " fragments"
puts "There are " + (all_frags_pos.flatten.length).to_s + " positions containing SNPs"

File.open("snp_pos.json", "w") do |f|
	f.write(all_frags_pos.to_json)
end

File.open("frags.json", "w") do |f|
	f.write(frags.to_json)
end

