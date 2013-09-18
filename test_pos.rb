#!/usr/bin/ruby

require "rubygems"
require "rinruby"

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
puts final_positions.length
puts unique_pos.length

