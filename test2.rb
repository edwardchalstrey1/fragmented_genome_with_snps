#!/usr/bin/ruby

require "rubygems"

require "json"

#a = ['x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x']
snp_pos = [8, 2, 11, 5, 0, 9, 3, 14, 15]
snp_pos = snp_pos.sort
f1 = ['x', 'x', 'x'] #can't touch this
f2 = ['x', 'x', 'x', 'x'] #can't touch this
f3 = ['x', 'x', 'x', 'x', 'x']
f4 = ['x', 'x', 'x', 'x']
d = [f1, f2, f3, f4]

d_ranges = []
d.each do |i|
	d_ranges << i.length + d_ranges[-1].to_i #adding the fragment lengths to get the upper bounds of ranges of positions on the original seq.
	# don't need to minus 1, because look below. we're doing < j
end
first_pos = []
first_pos << 0
first_pos << d_ranges[0..-2]
first_pos = first_pos.flatten
x = 0
p = 0
frags_pos = []
d_ranges.each do |j| #for each of the upper ranges (j) that the positions could be in
	frag_snp_pos = []
	# add all of the positions < j to a new array
	while snp_pos[x].to_i < j && !snp_pos[x].nil? do #make the loop quit before trying to access index of snp_pos that doesn't exist
		#puts "#{x} - #{snp_pos[x].to_i} - #{j}"
		frag_snp_pos << snp_pos[x].to_i 
		x += 1
	end
	z = 0
	y = first_pos[p].to_i 
	frag_snp_pos.each do |l|
		frag_snp_pos[z] = l - y #taking away the value at the first position of each fragments to get the true positions
		z += 1
	end
	p += 1
	frags_pos << frag_snp_pos #the position info is divided into the positions for each frag
end
puts snp_pos
puts
puts frags_pos
puts
puts first_pos


File.open("supertest.json", "w") do |f|
	f.write(frags_pos.to_json)
end


#puts "Same: " + (a == d.flatten).to_s


