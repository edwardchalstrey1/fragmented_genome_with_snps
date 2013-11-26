require "rubygems"
require "rinruby"

def density (snp_pos, length, n, m)
	length = (0..length-1).to_a
	ds2 = []
	y = 1
	length.each do |x| # the density values are inverted
		if snp_pos[y] != nil
			if x < snp_pos[0] # for bases before first SNP
				ds2 << n/(snp_pos[0]-x).to_f # take the base position from the first SNP position
			elsif x == snp_pos[0] # for the first SNP 
				ds2 << n/m # assign a value double that of the bases next to it
			elsif x < (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases less than halfway between the SNP y-1 and SNP y
				ds2 << n/(x-snp_pos[y-1]).to_f # take SNP y-1's position from the bases position
			elsif x < snp_pos[y] && x >= (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases more than halfway between the SNP y-1 and SNP y
				ds2 << n/(snp_pos[y]-x).to_f # take the bases position from SNP y's position
			elsif x == snp_pos[y] # for subsequent SNPs
				ds2 << n/m # assign a value double that of the bases next to it
			elsif x > snp_pos[y] # for bases higher than SNP y
				y+=1
				if snp_pos[y] != nil
					if x == snp_pos[y] # if the base is equal to the next SNP (y incremented)
						ds2 << n/m # assign a value double that of the bases next to it
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

def get_gradients (min_snps, n_div, m)
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
				y = density(snp_pos, length, n_div, m)
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

def example (frag_num, xs, ys, nu_ids, min_snps)
	min_snps_string =  min_snps.to_s
	z = nu_ids.index(frag_num)
	x = xs[z]
	y = ys[z]
	frag_num_string = frag_num.to_s
	myr = RinRuby.new(echo = false)
	myr.assign "x", x
	myr.assign "y", y
	myr.assign "frag_num_string", frag_num_string
	myr.assign "min_snps_string", min_snps_string
	myr.eval 'source("~/fragmented_genome_with_snps/skew_scatter.R")'
	myr.eval 'example(frag_num_string, x, y, min_snps_string)'
	myr.quit
end

def skew_scatta (min_snps, n_div, m, ex1, ex2, ex3) # n_div = number that is divided by the values of distance to get the 'density', m = value that n_div is divided by at SNP positions
	min_snps_string =  min_snps.to_s
	n_div_string = n_div.to_s
	gixs = get_gradients(min_snps, n_div, m)
	gradients = gixs[0]
	abs_gradients = gixs[4]
	nu_ids = gixs[1]
	xs = gixs[2]
	ys = gixs[3]
	fng = "grad_"+min_snps_string+"_"+n_div_string
	fna = "abs_"+min_snps_string+"_"+n_div_string
	scatta_dem(gradients, nu_ids, fng, "gradient", min_snps_string)
	scatta_dem(abs_gradients, nu_ids, fna, "absolute gradient", min_snps_string)
	example(ex1, xs, ys, nu_ids, min_snps)
	example(ex2, xs, ys, nu_ids, min_snps)
	example(ex3, xs, ys, nu_ids, min_snps)
end

Signal.trap("PIPE", "IGNORE")

#skew_scatta(2, 10000, 258, 681, 987)

skew_scatta(30, 10000, 0.5, 587, 694, 729)






