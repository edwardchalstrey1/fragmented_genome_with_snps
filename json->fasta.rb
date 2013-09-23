#!/usr/bin/ruby

require "rubygems"
require "json"
require 'fileutils'
require 'tempfile'

json = File.open('frags.json').read
frags = JSON.parse(json)

frag_strings = []
frags.each do |i|
	j = i.join.to_s #create a string of each sequence
	frag_strings << j #add these to a new array of the sequences in string format
end
frag_ids = []
x = 0
frag_strings.each do |f| #create an id for each of the frags
	frag_ids << ('>frag' + (x += 1).to_s)
end
frag_lengths = []
frags.each do |l|
	frag_lengths << l.length
end

json2 = File.open('snp_pos.json').read
snp_pos = JSON.parse(json2)

y = 0
snp_pos_strings = []
snp_pos.each do |n|
	s = "Length: " + (frag_lengths[y]).to_s + "  SNP number: " + n.length.to_s + "  SNP positions: " + n.to_s
	snp_pos_strings << s
	y += 1
end
id_and_pos = []
q = 0
snp_pos_strings.each do |h|
	one = frag_ids[q] + ", " + h
	id_and_pos << one
	q += 1
end
fastaformat_array = id_and_pos.zip(frag_strings)
fastaformat_array_shuf = fastaformat_array.shuffle
#puts fastaformat_array[0]
#puts fastaformat_array[1]
#puts
puts fastaformat_array_shuf[0]
#puts fastaformat_array_shuf[1]

File.open("frags.fasta", "w+") do |f|
	fastaformat_array_shuf.each { |i| f.puts(i) }
end
