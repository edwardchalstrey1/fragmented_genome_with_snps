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
		x = density(snp_pos, length, )
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

