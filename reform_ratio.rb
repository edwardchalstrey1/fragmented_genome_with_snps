require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'

def get_snp_data (vcf_file)
	vcfs_chrom = []
	vcfs_pos = []
	vcfs_info = []
	File.open(vcf_file, "r").each do |line| #get array of vcf lines, you can call a method on one line
		next if line =~ /^#/
		v = Bio::DB::Vcf.new(line)
		vcfs_chrom << v.chrom
		vcfs_pos << v.pos
		vcfs_info << v.info # so this will be an array of hashes of strings
	end
	num_snps_frag_hash = Hash.new(0)
	vcfs_chrom.each {|v| num_snps_frag_hash[v] +=1 } #we have the number of snps on each frag, by counting the repeats of each frag in the vcf
	#the frag_id(.chrom) is the key, the number of snps for that frag is the value. putting the number of snps for each frag into hash
	return vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info
end
def fasta_array (fasta_file)
	fasta = [] #we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
	Bio::FastaFormat.open(fasta_file).each do |i| #get array of fasta format frags, ##  WE NEED TO REORDER THE FASTA FRAGS HERE, TO TEST DIFFERENT ARRANGEMENTS
		fasta << i
	end
	return fasta
end
def snps_per_fasta_frag (snps_per_vcf_frag_hash, fasta_array)
	snps_per_frag_fasta_order = [] #use the id to identify the number of snps for that frag using the keys of snps_hash
	fasta_array.each do |frag|
		snps_per_frag_fasta_order << snps_per_vcf_frag_hash[frag.entry_id] #gives 0 for keys that don't exist = good, because the frags with 0 density would otherwise be missing
	end
	#now we have an array with the number of snps per frag in the same order as the fasta array, which we can get lengths from to calculate density
	return snps_per_frag_fasta_order
end
def get_positions (fasta, vcfs_chrom, vcfs_pos, snps_per_frag)
	pos = [] #get the snp positions for each frag, in an array of arrays
	n = 0
	fasta.each do |frag|
		x = 0
		each_fr_pos = []
		snps_per_frag[n].times do |j|
			if frag.entry_id == vcfs_chrom[x] #this assumes that frag_id == vcf.chrom then continues to for the number of snps (for that frag)
				each_fr_pos << vcfs_pos[x]
				x+=1
			else
				while frag.entry_id != vcfs_chrom[x]
					x+=1
				end
				each_fr_pos << vcfs_pos[x]
				x+=1
			end
		end
		pos << each_fr_pos #this gives empty arrays for frags with out snps, and a list of the positions of those with
		n+=1
	end
	return pos
end
def total_pos (pos, fasta_lengths) #pos is in the same order as the vcf, fasta lengths is in the fasta order, so this only works if they are in the correct order, 
	totals = []                    #but the totals are what we need for wrongly ordered frags too, they will be incorrect of course, but that is what we need to test
	x = 0						   # matching the ids for the info_hash and total_pos will work, because total pos is in the same order as pos, same order as VCF
	pos.each do |frag|
		if x == 0
			totals << frag.uniq
			x+=1
		else
			tot_frag = []
			lengths = []
			fasta_lengths[0..x-1].each do |p|
				lengths << p
			end
			so_far = lengths.inject(:+) # this needs to be the length of the frags, not the number of snps
			frag.each do |i|
				tot_frag << so_far-1 + i
			end
			totals << tot_frag.uniq
			x+=1
		end
	end
	return totals
end
def het_hom (actual_pos, vcfs_info) #actual_pos in same order as VCF
	het = []
	hom = []
	x = 0
	actual_pos.flatten.each do |snp|
		if vcfs_info[x] == {"AF"=>"1.0"} # homozygous SNPs have AF= 1.0, we can change this to a range for real data
			hom << snp
		elsif vcfs_info[x] == {"AF"=>"0.5"}
			het << snp
		end
		x+=1
	end
	return het, hom
end
def prime? (n)
	for d in 2..(n - 1)
		if (n % d) == 0
    		return false
    	end
	end
	true
end
def division (frags, prime) #number of frags
	x = 1.5
	if prime == "n"
		until frags.length/x == (frags.length/x).to_i && x == x.to_i && x <= frags.length
			x = (frags.length/10).to_f + rand(frags.length).to_f
		end
	elsif prime == "p"
		until frags.length/x == (frags.length/x).to_i && x == x.to_i && x < frags.length
			x = (frags.length/10).to_f + rand(frags.length).to_f
		end
	end		
	return x
end
def recombine (mum, dad)
	kid = []
	1.times do
		mum1 = "bogey"
		x = division(mum, "n")
		if x == mum.length && prime?(x) == false
			#puts "redoing1..."
			redo
		elsif x == mum.length # to compensate for datasets with a prime number of fragments:
			#puts "we gotta prime"
			ig = rand(mum.length)-1 # choose a random element of the fasta array to ignore
			mig = mum[ig] # we can add these frags back at their original positions after recombination
			dig = dad[ig]
			mum1 = mum.dup
			dad1 = dad.dup
			mum1.delete_at(ig)
			dad1.delete_at(ig)
			x = division(mum1, "p")
		end
		if mum1 != "bogey"
			mum1 = mum1.each_slice(x).to_a
			dad1 = dad1.each_slice(x).to_a
		else
			mum1 = mum.each_slice(x).to_a
			dad1 = dad.each_slice(x).to_a
		end
		#puts "Ignored a frag? " +(ig != nil).to_s
		# Let's say we use one chunk of the dad solution, the rest mum
		#puts "x is " +x.to_s
		ch = rand(dad1.length)-1 # choose one of the chunks of fragments to keep from dad
		child = mum1.dup.flatten
		y = 0
		pos_array = []
		mum1[ch].each do |frag| # place each frag in the equivalent mum chunk into the position it's corresponding frag (from dad) occupies in mum
			chunk_frag = dad1[ch][y] # the equivalent frag in the chosen dad chunk
			pos = mum1.flatten.index(chunk_frag) # the position of the dad chunk frag in mum
			c_pos = mum1.flatten.index(frag) # the position of the frag in mum
			pos_array << pos
			#puts pos
			y+=1
		end
		if pos_array.include?(nil)
			#puts "redoing2..."
			redo
		else
			y = 0
			pos_array.each do |pos|
				unless dad1[ch].include?(mum1[ch][y])
					child[pos] = mum1[ch][y]
					child[mum1.flatten.index(mum1[ch][y])] = dad1[ch][y] # swapping the positions of the frag and chunk frag, to give their positions in child
				end
				y+=1
			end
		end
		#puts "kid length same: "+(mum1.flatten.length == child.length).to_s
		if ig != nil
			if dad1[ch].include?(dig) # add the ignored fragment at it's 
				child.insert(ig, dig)
			else
				child.insert(ig, mig)
			end
		end
		#puts "kid unique? "+(child == child.uniq).to_s
		kid << child
	end
	return kid[0]
end
def fasta_child (fasta)
	mum = fasta.shuffle
	dad = fasta.shuffle
	child = recombine(mum, dad)
	return child
end
def fitness (fasta, vcf)
	fasta_ids = []
	fasta_lengths = []
	fasta.each do |i|
		fasta_ids << i.entry_id
		fasta_lengths << i.length
	end
	snp_data = get_snp_data(vcf)
	vcfs_chrom = snp_data[0] #array of vcf frag ids
	vcfs_pos = snp_data[1] #array of all the snp positions (fragments with snps)
	snps_hash = snp_data[2] #hash of each fragment from vcf, and it's number of snps
	vcfs_info = snp_data[3]
	snps_per_frag = snps_per_fasta_frag(snps_hash, fasta) #array of no. of snps per frag in same order as fasta
	pos = get_positions(fasta, vcfs_chrom, vcfs_pos, snps_per_frag) #get snp positions for each frag in array of arrays
	actual_pos = total_pos(pos, fasta_lengths)
	het_hom_snps = het_hom(actual_pos, vcfs_info)
	het = het_hom_snps[0]
	hom = het_hom_snps[1]
	myr = RinRuby.new(echo=false)
	myr.assign "het_snps", het
	myr.assign "hom_snps", hom
	myr.eval "source('~/fragmented_genome_with_snps/ratio.R')"
	coeff = myr.pull "cor(qqp$x,qqp$y)"
	myr.quit
	#write_txt("arabidopsis_datasets/"+ARGV[0].to_s+"/het_snps", het)
	#write_txt("arabidopsis_datasets/"+ARGV[0].to_s+"/hom_snps", hom)
	return coeff
end

### Genetic Algorithm

def initial_population(fasta, size)
	population = []
	size.times do
		chromosome = fasta.shuffle
		population << chromosome
	end
	return population
end
def select(pop, vcf)
	fits = []
	pop.each do |sol| #solution
		fits << fitness(sol, vcf)
	end
	pop_fits = fits.zip(pop).sort
	return pop_fits[-1], pop_fits[-2]
end
def new_population(parent1, parent2, size)
	pop = []
	size-2.times do
		pop << recombine(parent1, parent2)
	end
	return pop
end
def evolve(fasta, vcf, gen, pop_size)
	pop = initial_population(fasta, pop_size)
	best2 = select(pop, vcf)
	puts
	puts "Gen 0"
	puts "Coefficient 1best= "+(best2[0][0]).to_s
	puts "Coefficient 2best= "+(best2[1][0]).to_s
	puts
	y=1
	gen.times do
		pop = new_population(best2[0][1], best2[1][1], pop_size)
		best2 = select(pop, vcf)
		puts "Gen"+y.to_s
		puts "Coefficient 1best= "+(best2[0][0]).to_s
		puts "Coefficient 2best= "+(best2[1][0]).to_s
		puts
		y+=1
	end
end

vcf = 'arabidopsis_datasets/'+ARGV[0].to_s+'/snps.vcf'
fasta = fasta_array('arabidopsis_datasets/'+ARGV[0].to_s+'/frags_shuffled.fasta') #array of fasta format fragments, and entry_ids

evolve(fasta, vcf, 10, 20)
