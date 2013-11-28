require "rubygems"
require "rinruby"

def density (snp_pos, length, n, m)
	length = (0..length-1).to_a
	ds = []
	if snp_pos.length != 1
		y = 1
		length.each do |x| # the density values are inverted
			if snp_pos[y] != nil
				if x < snp_pos[0] # for bases before first SNP
					ds << n/(snp_pos[0]-x).to_f # take the base position from the first SNP position
				elsif x == snp_pos[0] # for the first SNP 
					ds << n/m # assign a value double that of the bases next to it
				elsif x < (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases less than halfway between the SNP y-1 and SNP y
					ds << n/(x-snp_pos[y-1]).to_f # take SNP y-1's position from the bases position
				elsif x < snp_pos[y] && x >= (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases more than halfway between the SNP y-1 and SNP y
					ds << n/(snp_pos[y]-x).to_f # take the bases position from SNP y's position
				elsif x == snp_pos[y] # for subsequent SNPs
					ds << n/m # assign a value double that of the bases next to it
				elsif x > snp_pos[y] # for bases higher than SNP y
					y+=1
					if snp_pos[y] != nil
						if x == snp_pos[y] # if the base is equal to the next SNP (y incremented)
							ds << n/m # assign a value double that of the bases next to it
						else
							ds << n/(x-snp_pos[y-1]).to_f # if not, take the SNP y-1's position from the bases position
						end
					else
						ds << n/(x-snp_pos[y-1]).to_f # if there are no SNP's left, take the SNP y-1's position from the bases position
					end
				end
			else
				ds << n/(x-snp_pos[-1]).to_f
			end
		end
	else
		length.each do |x|
			if x < snp_pos[0]
				ds << n/(snp_pos[0]-x)
			elsif x == snp_pos[0]
				ds << n/m
			elsif x > snp_pos[0]
				ds << n/(x-snp_pos[0])
			end
		end
	end
	return ds
end

def density2 (snp_pos, length, n)
	length = (0..length-1).to_a
	ds = []
	if snp_pos.length != 1
		y = 1
		length.each do |x| # the density values are inverted
			if snp_pos[y] != nil
				if x < snp_pos[0] # for bases before first SNP
					ds << n - (snp_pos[0]-x).to_f # take the base position from the first SNP position
				elsif x == snp_pos[0] # for the first SNP 
					ds << n # assign a high value
				elsif x < (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases less than halfway between the SNP y-1 and SNP y
					ds << n - (x-snp_pos[y-1]).to_f # take SNP y-1's position from the bases position
				elsif x < snp_pos[y] && x >= (snp_pos[y]+snp_pos[y-1]).to_f/2 # for bases more than halfway between the SNP y-1 and SNP y
					ds << n - (snp_pos[y]-x).to_f # take the bases position from SNP y's position
				elsif x == snp_pos[y] # for subsequent SNPs
					ds << n # assign a high value
				elsif x > snp_pos[y] # for bases higher than SNP y
					y+=1
					if snp_pos[y] != nil
						if x == snp_pos[y] # if the base is equal to the next SNP (y incremented)
							ds << n # assign a high value
						else
							ds << n - (x-snp_pos[y-1]).to_f # if not, take the SNP y-1's position from the bases position
						end
					else
						ds << n - (x-snp_pos[y-1]).to_f # if there are no SNP's left, take the SNP y-1's position from the bases position
					end
				end
			else
				ds << n - (x-snp_pos[-1]).to_f
			end
		end
	else
		length.each do |x|
			if x < snp_pos[0]
				ds << n - (snp_pos[0]-x)
			elsif x == snp_pos[0]
				ds << n
			elsif x > snp_pos[0]
				ds << n - (x-snp_pos[0])
			end
		end
	end
	return ds
end

def get_gradients (min_snps, n_div, m, d)
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
				if d == 1
					y = density(snp_pos, length, n_div, m)
				elsif d == 2
					y = density2(snp_pos, length, n_div)
				end	
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

def example (frag_num, xs, ys, nu_ids, min_snps, d, m_str)
	min_snps_string =  min_snps.to_s
	z = nu_ids.index(frag_num)
	x = xs[z]
	y = ys[z]
	frag_num_string = frag_num.to_s
	d_str = 'd'+d.to_s
	myr = RinRuby.new(echo = false)
	myr.assign "m", m_str
	myr.assign "d", d_str
	myr.assign "x", x
	myr.assign "y", y
	myr.assign "frag_num_string", frag_num_string
	myr.assign "min_snps_string", min_snps_string
	myr.eval 'source("~/fragmented_genome_with_snps/skew_scatter.R")'
	myr.eval 'example(frag_num_string, x, y, min_snps_string, d, m)'
	myr.quit
end

def skew_scatta (min_snps, n_div, m, d, ex1, ex2, ex3) 
	min_snps_string =  min_snps.to_s
	n_div_string = n_div.to_s
	m_str = m.to_s
	gixs = get_gradients(min_snps, n_div, m, d)
	gradients = gixs[0]
	abs_gradients = gixs[4]
	nu_ids = gixs[1]
	xs = gixs[2]
	ys = gixs[3]
	if d == 1
		fng = "grad_"+min_snps_string+"_"+n_div_string+"_d1_m"+m_str
		fna = "abs_"+min_snps_string+"_"+n_div_string+"_d1_m"+m_str
	elsif d == 2
		fng = "grad_"+min_snps_string+"_"+n_div_string+"_d2"
		fna = "abs_"+min_snps_string+"_"+n_div_string+"_d2"
	end
	scatta_dem(gradients, nu_ids, fng, "gradient", min_snps_string)
	scatta_dem(abs_gradients, nu_ids, fna, "absolute gradient", min_snps_string)
	example(ex1, xs, ys, nu_ids, min_snps, d, m_str)
	example(ex2, xs, ys, nu_ids, min_snps, d, m_str)
	example(ex3, xs, ys, nu_ids, min_snps, d, m_str)
end
# n_div = number that is divided by the values of distance to get the 'density',
# or the number that the distance is taken from in density2.
# m = value that n_div is divided by at SNP positions
# d = density method
Signal.trap("PIPE", "IGNORE")

#skew_scatta(2, 10000, 0.5, 1, 258, 681, 987)
#skew_scatta(30, 10000, 0.5, 1, 587, 694, 729)
#skew_scatta(1, 10000, 0.5, 1, 258, 681, 987)
#skew_scatta(1, 20000, 0.5, 2, 258, 681, 987)
#skew_scatta(30, 20000, 0.5, 2, 587, 694, 729)
#skew_scatta(1, 10000, 0.25, 1, 258, 681, 729)


def how_scatta (id) # id = frag no.
	lengths = []
	File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_fasta_lengths.txt").each {|line| lengths << line.to_i}
	ids = []
	File.open("arabidopsis_datasets/dataset5/skew_scatter/ex_ids_w_snps.txt").each {|line| ids << line.to_i}
	z = ids.index(id)
	length = lengths[z]
	snp_pos = []
	File.open("arabidopsis_datasets/dataset5/skew_scatter/snps"+id.to_s+".txt").each {|line| snp_pos << line.to_i}
	y = []
	snp_pos.length.times {|i| (y << 1)}
	myr = RinRuby.new(echo=false)
	myr.assign 'length', length
	myr.assign 'y', y
	myr.assign 'frag_num_string', id.to_s
	myr.assign 'snp_pos', snp_pos
	myr.eval 'source("~/fragmented_genome_with_snps/skew_scatter.R")'
	myr.eval 'how_scatta(frag_num_string, snp_pos, y, length)'
	myr.quit
end

#1234 frags
#how_scatta(335)
#how_scatta(658)
#how_scatta(659)
#how_scatta(707)
#how_scatta(842)
#how_scatta(1019)









