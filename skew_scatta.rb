require "rubygems"
require "rinruby"

lengths = []
File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_fasta_lengths.txt").each {|line| lengths << line}
ids = []
File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_ids_w_snps.txt").each {|line| ids << line}
gradients = []
xs = []
ys = []
q = 1
n = 1
nu_ids = []
lengths.each do |length|
	snp_pos = []
	File.open("arabidopsis_datasets/dataset5/skew_scatter/snps"+ids[q].to_s+".txt").each {|line| snp_pos << line}
	if snp_pos.length > 6
		nu_ids << ids[q]
		#make x and y
		xs << x
		ys << y
		n+=1
		myr = RinRuby.new(echo = false)
		myr.assign "x", x
		myr.assign "y", y
		myr.eval "gradient <- (coef(lm(y ~ x)))[2]"
		gradients << (myr.pull "gradient")
	end
	q+=1
end




length = 20
snp_pos = [5,11,14,16,17,19]
ds = []
ds << 10/(snp_pos[1]-snp_pos[0]).to_f
n = 2
m = 0
snp_pos[1..-2].each do |pos|
	ds << 10/(((pos-snp_pos[m])+(snp_pos[n]-pos)).to_f/2)
	n+=1
	m+=1
end
ds << 10/(snp_pos[-1]-snp_pos[-2]).to_f
#puts ds



y = 1
ds2 = []
length = (1..length).to_a
length.each do |x|
	if x < snp_pos[0]
		ds2 << (snp_pos[y]-x).to_f
	elsif x == snp_pos[0]
		ds2 << 1.to_f
	elsif x < snp_pos[y]
		ds2 << ((x-snp_pos[y-1])+(snp_pos[y]-x)).to_f/2
	elsif x == snp_pos[y]
		ds2 << 1
	else
		y+=1
		ds2 << ((x-snp_pos[y-1])+(snp_pos[y]-x)).to_f/2
	end
end
puts ds2
