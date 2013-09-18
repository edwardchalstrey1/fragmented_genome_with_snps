#!/usr/bin/ruby

require "rubygems"
require "rinruby"

#Create array of SNP positions:
positions = []
while positions.length < 1 do # I DON'T UNDERSTAND WHY THIS WORKS, BUT IT DOES
	myr = RinRuby.new(echo = false) #new instance of rinruby, with supressed output from R
	myr.eval "x <- rpois(1000, 5000000)" #poisson distibution of 1000 snp positions about a mean of 5m
	snp_pos = myr.pull "x" #pulling the position data to an array
	snp_pos.each{
		|i|
	
		if positions.include? i
			next
		else
			positions << i
		end
	}
end #^adding the unique positions (935/1000) to positions[]

#Create sequence with SNPs at designated positions, and a reference sequence with no SNPs
a = ['a', 'g', 't', 'g', 'a']
ref = a*2000000 # reference sequence one base per array slot
snp_seq = a*2000000
positions.each{
	|i|
	snp_seq[i] = 'c' #replacing the slots at the snp positions with snp's (always c for simplicity)
}

#THE LOOP BELOW WORKS, BUT IS VERY VERY SLOW
#Split sequences into random sized fragments
class Fragment

	def initialize(seq)
		fragment(seq)
	end

	def fragment(seq)
		frags = []
		until frags.flatten.join.to_s.split(//).length >= 10000000 do #once the frags array has all 50 bases from snp_seq
			frags << (seq.each_slice(rand(2000000) + 500000).to_a)[0] #give one frag of 2-5b to the frags array
			seq_string = seq.join.to_s  #create a string from the sequence for slicing
			frags_string = frags[-1].join.to_s #create a string of the most recently added fragment to slice off
			seq_string.slice! frags_string #slice the fragment from the start of the string and... 
			seq = seq_string.split(//) #modify the sequence so it no longer contains the removed frags, then the loop can repeat
			puts frags.flatten.join.to_s.split(//).length
		end
		puts frags.length
	end
end

ref_frags = Fragment.new(ref)
snp_frags = Fragment.new(snp_seq)
