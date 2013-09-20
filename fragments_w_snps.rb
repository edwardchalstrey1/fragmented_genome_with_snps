#!/usr/bin/ruby

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
#puts snp_counts#.length

snp_pos_length = []
snp_pos = []
x = Random.new
snp_counts.each do |l| #for each count do
	snp_breaks.each do |i| # for each break: split the breaks to get the ranges i..i+1
		ii = i.to_i
		iii = (i.to_i) + 9999 #THIS BIT IS STILL KINDA CHEATING
		snpo = []
		(ii..iii).each do |j| #for each range WHY ARE WE DOING 10K EACH TIME HERE?: we need to do a random no. for each position
			newr = x.rand(1..200000) #make a random no. within range 
			if newr <= l
				snpo << j #add the position (j) to array if its random no. is <= count
			end 
		end #snpo is now an array of snp positions for the range i to i+1
		snp_pos << snpo #snp_pos is now an array snpo arrays for each i
	end
	snp_pos_length << snp_pos.length # because of loop, length increases, making array of each new length for each iteration
end #snp_pos now contains (the snpo arrays for each i)*l 
length_index = (1..snp_pos_length.length).to_a
a = [length_index, snp_pos_length] 
snpo_positions = a.transpose.map {|x| x.reduce(:+)} 
#adding the snp_pos lengths to their indices(+1), to find the relevant snpo's. for each count(l), 
#we need the lth snpo in snp_pos, as this contains the positions for the break range associated with that count:
snpo_positions.delete_at(-1) #last element needs to be removed, all the others contain the snpo position of relevance
actual_pos = []
actual_pos << snp_pos[0] #adding the 0th snpo, the snp positions for l=0 i=0.
snpo_positions.each do |c| # then adding the rest, so only the ones at the specified positions: l=x i=x
	actual_pos << snp_pos[c]
end
all_snp_positions = actual_pos.flatten
#puts all_snp_positions.length
myr.assign "pos_l", all_snp_positions
myr.eval "pos_s <- sample(pos_l, 1000)"
final_positions = myr.pull "pos_s"

unique_pos = []
final_positions.each{
	|q|
	
	if unique_pos.include? q
		next
	else
		unique_pos << q
	end
}
puts "There are " + final_positions.length.to_s + " SNP\'s"
puts unique_pos.length.to_s + " of the SNP positions are unique" #this proves each snp position is unique
puts
#Create sequence with SNPs at designated positions
a = ['a', 'g', 't', 'g', 'a']
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
puts "The sequence array == the flattened fragment array: " + (snp_seq == frags.flatten).to_s

#Find the positions of snps on the fragments, based on the positions in the un-fragmented sequence


File.new("frags.json", "w")
File.open("frags.json", "w") do |f|
	f.write(frags.to_json)
end

