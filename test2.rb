def density (snp_pos, length, n)
	length = (1..length).to_a
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
		end
	end
	return ds2
end


length = 20
snp_pos = [5,11,14,16,17]

puts density(snp_pos, length, 10)