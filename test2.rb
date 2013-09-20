a = ['x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x', 'x']
snp_pos = [0, 2, 3, 5, 8, 9, 11]
f1 = ['x', 'x', 'x'] #can't touch this
f2 = ['x', 'x', 'x', 'x'] #can't touch this
f3 = ['x', 'x', 'x', 'x', 'x']
d = [f1, f2, f3]

d_ranges = []
d.each do |i|
	d_ranges << i.length + d_ranges[-1].to_i
end
x = 0
#puts snp_pos[x]
frags_pos = []
d_ranges.each do |j|
	frag_snp_pos = []
	while snp_pos[x].to_i < j && !snp_pos[x].nil? do #make the loop quit before trying to access index of snp_pos that doesn't exist
		puts "#{x} - #{snp_pos[x].to_i} - #{j}"
		frag_snp_pos << snp_pos[x].to_i
		x += 1
		puts x
	end
	frags_pos << frag_snp_pos
end




#puts "Same: " + (a == d.flatten).to_s


