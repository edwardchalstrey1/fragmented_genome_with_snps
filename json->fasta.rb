#!/usr/bin/ruby

require "rubygems"
require "json"
require 'fileutils'
require 'tempfile'

json = File.open('frags.json').read
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
#puts fastaformat_array

File.new("frags.fasta", "w+")
File.open("frags.fasta", "w+") do |f|
	fastaformat_array.each { |i| f.puts(i) }
end
