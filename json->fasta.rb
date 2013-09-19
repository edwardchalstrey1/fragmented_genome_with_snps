#!/usr/bin/ruby

require "rubygems"
require "json"
require 'fileutils'
require 'tempfile'

json = File.open('test.json').read
frags = JSON.parse(json)

frag_strings = []
frags.each do |i|
	j = i.join.to_s
	frag_strings << j
end
#puts frag_strings

frag_ids = []
x = 0
frag_strings.each do |f|
	frag_ids << ('>frag' + (x += 1).to_s)
end
fastaformat_array = frag_ids.zip(frag_strings).flatten.compact
puts fastaformat_array

y = -1
File.new("test.fasta", "w+")
File.open("test.fasta", "w+") do |i|
	i.write(fastaformat_array)
end
