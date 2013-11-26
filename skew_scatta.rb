require "rubygems"
require "rinruby"

def density (snp_pos, length, n)
	length = (0..length-1).to_a
	ds2 = []
	y = 1
	length.each do |x| # the density values are inverted
		if snp_pos[y] != nil
			if x < snp_pos[0] # for bases before first SNP
				ds2 << n/(snp_pos[0]-x).to_f # take the base position from the first SNP position
			elsif x == snp_pos[0] # for the first SNP 
				ds2 << n/0.5 # assign a value double that of the bases next to it
			elsif x < (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases less than halfway between the SNP y-1 and SNP y
				ds2 << n/(x-snp_pos[y-1]).to_f # take SNP y-1's position from the bases position
			elsif x < snp_pos[y] && x >= (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases more than halfway between the SNP y-1 and SNP y
				ds2 << n/(snp_pos[y]-x).to_f # take the bases position from SNP y's position
			elsif x == snp_pos[y] # for subsequent SNPs
				ds2 << n/0.5 # assign a value double that of the bases next to it
			elsif x > snp_pos[y] # for bases higher than SNP y
				y+=1
				if snp_pos[y] != nil
					if x == snp_pos[y] # if the base is equal to the next SNP (y incremented)
						ds2 << n/0.5 # assign a value double that of the bases next to it
					else
						ds2 << n/(x-snp_pos[y-1]).to_f # if not, take the SNP y-1's position from the bases position
					end
				else
					ds2 << n/(x-snp_pos[y-1]).to_f # if there are no SNP's left, take the SNP y-1's position from the bases position
				end
			end
		else
			ds2 << n/(x-snp_pos[-1]).to_f
		end
	end
	return ds2
end

def get_gradients (min_snps, n_div)
	lengths = []
	File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_fasta_lengths.txt").each {|line| lengths << line.to_i}
	ids = []
	File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_ids_w_snps.txt").each {|line| ids << line.to_i}
	gradients = []
	xs = []
	ys = []
	q = 1
	n = 1
	nu_ids = []
	lengths.each do |length|
		snp_pos = []
		if ids[q] != nil
			File.open("arabidopsis_datasets/dataset5/skew_scatter/snps"+ids[q].to_s+".txt").each {|line| snp_pos << line.to_i}
			if snp_pos.length >= min_snps
				nu_ids << ids[q]
				y = density(snp_pos, length, n_div)
				x = (0..length-1).to_a
				xs << x
				ys << y
				n+=1
				myr = RinRuby.new(echo = false)
				myr.assign "x", x
				myr.assign "y", y
				myr.eval "gradient <- (coef(lm(y ~ x)))[2]"
				gradient = (myr.pull "gradient")
				gradients << gradient
				myr.quit
			end
			puts q
			q+=1
		end
	end
	abs_gradients = []
	gradients.each {|g| abs_gradients << g.abs}
	return gradients, nu_ids, xs, ys, abs_gradients
end

def scatta_dem (gradients, nu_ids, filename, grad_string, min_snps_string)
	myr = RinRuby.new(echo = false)
	myr.assign "gradients", gradients
	myr.assign "nu_ids", nu_ids
	myr.assign "filename", filename
	myr.assign "grad_string", grad_string
	myr.assign "min_snps_string", min_snps_string
	myr.eval 'source("~/fragmented_genome_with_snps/skew_scatter.R")'
	myr.eval 'scatta(gradients, nu_ids, filename, grad_string, min_snps_string)'
	myr.quit
end

def example (frag_num, xs, ys)
	x = xs[frag_num-1]
	y = ys[frag_num-1]
	frag_num_string = frag_num.to_s
	myr = RinRuby.new(echo = false)
	myr.assign "x", x
	myr.assign "y", y
	myr.assign "frag_num_string", frag_num_string
	myr.eval 'source("~/fragmented_genome_with_snps/skew_scatter.R")'
	myr.eval 'example(frag_num_string, x, y)'
	myr.quit
end
Signal.trap("PIPE", "IGNORE")
gixs = get_gradients(2, 10000)
gradients = gixs[0]
abs_gradients = gixs[4]
nu_ids = gixs[1]
xs = gixs[2]
ys = gixs[3]
scatta_dem(gradients, nu_ids, "grad_2_10K", "gradient", "2")
scatta_dem(abs_gradients, nu_ids, "abs_2_10K", "absolute gradient", "2")
example(258, xs, ys)
example(681, xs, ys)
example(974, xs, ys)

#snp_pos = []
#File.open("arabidopsis_datasets/dataset5/skew_scatter/snps230.txt").each {|line| snp_pos << line.to_i}
#lengths = []
#File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_fasta_lengths.txt").each {|line| lengths << line.to_i}
#l = lengths[229]
#d = (density(snp_pos, l, 10000))
#myr = RinRuby.new(echo = false)
#myr.assign "d", d
#myr.eval "s <- length(d)"
#s = myr.pull "s"
#puts s








