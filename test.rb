#!/usr/bin/ruby

require "rubygems"
require "json"

a = [1,2,3,4,5,6,7,8,9,10]
b = ['a','b','c', 'g', 'a', 'a', 'a', 'k']
c = [a, b]
all = c.flatten
x = all.length - 17
(c[-1])[-(x)..-1] = 'p'
c[-1].delete_if {|i| i > 'o'}

puts c

File.new("test.json", "w")

File.open("test.json", "w") do |f|
	f.write(c.to_json)
end