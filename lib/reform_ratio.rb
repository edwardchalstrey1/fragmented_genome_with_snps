#require 'rubygems'
#require 'bio-samtools'
#require 'bio'
#require 'rinruby' ## These are needed in a script that uses this class

class ReformRatio

# Input: Array of Bio::FastaFormat entries
# Output 0: Array of identifiers
# Output 1: Array of lengths (integers)
def self.fasta_id_n_lengths (fasta)
	ids = []
	lengths = []
	fasta.each do |i|
		ids << i.entry_id
		lengths << i.length
	end
	return ids, lengths
end

# Input: VCF file
# Output 0: Array of VCF chrom field (fragment identifiers)
# Output 1: Array of VCF pos field (positions of snps on fragments)
# Output 2: Hash with fragment id keys and the corresponding number of snps as an integer values
# Output 3: Array of VCF info field (hashes of the info values e.g. key: AF, value: allele frequency)
def self.get_snp_data (vcf_file)
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

# Input: FASTA file
# Output: Array of Bio::FastaFormat entries
def self.fasta_array (fasta_file)
	fasta = [] #we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
	Bio::FastaFormat.open(fasta_file).each do |i| #get array of fasta format frags, ##  WE NEED TO REORDER THE FASTA FRAGS HERE, TO TEST DIFFERENT ARRANGEMENTS
		fasta << i
	end
	return fasta
end

# Input 0: Hash with fragment id keys and the corresponding number of snps as an integer values
# Input 1: Array of Bio::FastaFormat entries
# Output: Array of the number of snps per fragment, in the same order as the input 1 fasta array
def self.snps_per_fasta_frag (snps_per_vcf_frag_hash, fasta_array)
	snps_per_frag_fasta_order = [] #use the id to identify the number of snps for that frag using the keys of snps_hash
	fasta_array.each do |frag|
		snps_per_frag_fasta_order << snps_per_vcf_frag_hash[frag.entry_id] #gives 0 for keys that don't exist = good, because the frags with 0 density would otherwise be missing
	end
	#now we have an array with the number of snps per frag in the same order as the fasta array
	return snps_per_frag_fasta_order
end

# Input 0: Array of Bio::FastaFormat entries
# Input 1: Array of VCF chrom field (fragment identifiers)
# Input 2: Array of VCF pos field (positions of snps on fragments)
# Input 3: Array of the number of snps per fragment, in the same order as the input 0 fasta array
# Input 4: Array of VCF info field (hashes of the info values e.g. key: AF, value: allele frequency, single key/value)
# Output 0: The snp positions for each frag, in an array of arrays (each sub array contains the snp positions for one frag, and the sub arrays are ordered according to the order of the input 1 fasta)
# Output 1: The info hashes (of each snp, single key/value) for each frag, in an array of arrays (each sub array contains the info hashes for one frag, and the sub arrays are ordered according to the order of the input 1 fasta)
def self.get_positions (fasta, vcfs_chrom, vcfs_pos, snps_per_frag, vcfs_info)
	pos = [] #get the snp positions for each frag, in an array of arrays
	info = []
	n = 0
	fasta.each do |frag|
		x = 0
		each_fr_pos = []
		each_fr_info = []
		snps_per_frag[n].times do |j|
			if frag.entry_id == vcfs_chrom[x] #this assumes that frag_id == vcf.chrom then continues to for the number of snps (for that frag)
				each_fr_pos << vcfs_pos[x]
				each_fr_info << vcfs_info[x]
				x+=1
			else
				while frag.entry_id != vcfs_chrom[x]
					x+=1
				end
				each_fr_pos << vcfs_pos[x]
				each_fr_info << vcfs_info[x]
				x+=1
			end
		end
		pos << each_fr_pos #this gives empty arrays for frags with out snps, and a list of the positions of those with
		info << each_fr_info
		n+=1
	end
	return pos, info
end

# Input 0: The snp positions for each frag, in an array of arrays (for a given fragment permutation)
# Input 1: Array of fragment lengths (integers) (for the same fragment permutation as Input 0)
# Output: Array of the snp positions in the genome (assuming the genome is ordered according to the fragment permutation in Input 0)
def self.total_pos (pos, fasta_lengths)  # both args in same order as fasta = good. this all works!
	totals = []                    
	x = 0						   
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
	return totals.flatten
end

# Input 0: Array of the snp positions in the genome
# Input 1: Array of VCF info field (hashes of the info values e.g. key: AF, value: allele frequency)
# Output 0: Array of all the heterozygous snp positions
# Output 1: Array of all the homozygous snp positions
def self.het_hom (actual_pos, vcfs_info) #actual_pos in same order as fasta perm. now so is info
	het = []
	hom = []
	x = 0
	actual_pos.each do |snp|
		if vcfs_info.flatten[x] == {"AF"=>"1.0"} # homozygous SNPs have AF= 1.0, we can change this to a range for real data
			hom << snp
		elsif vcfs_info.flatten[x] == {"AF"=>"0.5"}
			het << snp
		end
		x+=1
	end
	return het, hom
end

# Input 0: Array of fragment identifiers in the correct order
# Input 1: Array of fragment identifiers in the order you wish to check
# Output 0: Ordinal similarity score value (0 = correct order)
def self.rearrangement_score (frags_original_order, rearranged)
	position_each_frag_id_in_d = frags_original_order.map{|x| rearranged.index(x)} #works out the index of frags_original_order values in rearranged
	index_values = Array(0..(frags_original_order.length - 1)) # index values that frags_original_order originally at
	both = []
	both << position_each_frag_id_in_d
	both << index_values
	difference = both.transpose.map {|x| x.reduce(:-)} # taking away old position from new position, to find the distance that the frag has moved when re-ordered
	difference_abs = []
	difference.each do |i|
		difference_abs << i.abs
	end
	score = difference_abs.inject(:+) #high score = bad, score of 0 means the fragments in the right order
	return score
end

# Input 0: Filename by which to save an array to a .txt file, one value per line
# Input 1: Array to save
def self.write_txt (filename, array)
	File.open(filename, "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end

# Input: An integer you wish to know whether it's a prime
# Output: true/false
def self.prime? (n)
	for d in 2..(n - 1)
		if (n % d) == 0
    		return false
    	end
	end
	true
end

# Genetic Algorithm
#------------------

# Input: Array
# Output: A random integer that the length of the Input 0 array can be divided by to get another integer (the randomly chosen size of chunks that permutations will be split into, in the recombine/mutate methods)
def self.division (frags) #number of frags must be > 10
	x = 1.5
	until frags.length/x == (frags.length/x).to_i && x == x.to_i && x <= frags.length
		x = (frags.length/10).to_f + rand(frags.length).to_f
	end
	return x
end

# Input 0: A parent permutation array of Bio::FastaFormat entries (or any array of unique objects)
# Input 1: A second parent permutation array of Bio::FastaFormat entries (or any array of the same unique objects as input 0)
# Output: A child permutation array of Bio::FastaFormat entries, whose order is a recombination of the parent permutations (the same unique objects ordered differently to either input)
def self.recombine (mum, dad)
	kid = []
	1.times do # so we can use redo
		mum1 = "bogey"
		x = division(mum)
		if x == mum.length && prime?(x) == false
			redo
		elsif x == mum.length # to compensate for datasets with a prime number of fragments:
			ig = rand(mum.length)-1 # choose a random element of the fasta array to ignore
			mig = mum[ig] # we can add these frags back at their original positions after recombination
			dig = dad[ig]
			mum1 = mum.dup
			dad1 = dad.dup
			mum1.delete_at(ig)
			dad1.delete_at(ig)
			x = division(mum1)
		end
		if mum1 != "bogey"
			mum1 = mum1.each_slice(x).to_a
			dad1 = dad1.each_slice(x).to_a
		else
			mum1 = mum.each_slice(x).to_a
			dad1 = dad.each_slice(x).to_a
		end # Let's say we use one chunk of the dad solution, the rest mum
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
		if ig != nil
			if dad1[ch].include?(dig) # add the ignored fragment at it's 
				child.insert(ig, dig)
			else
				child.insert(ig, mig)
			end
		end
		if child == mum || child == dad || child != child.uniq
			redo
		end
		kid << child #so we can access this outside the loop
	end
	return kid[0]
end

# Input: A permutation array of Bio::FastaFormat entries (or any array)
# Output: The input permutation array of Bio::FastaFormat entries, with a small change in the fragment order
def self.mutate (fasta)
	x = 0
	until x > 2
		x = division(fasta)
	end
	sliced = fasta.each_slice(x).to_a
	e = rand(sliced.length-1).to_i
	sliced[e] = sliced[e].shuffle
	return sliced.flatten
end

# Input: A permutation array of Bio::FastaFormat entries
# Output: The input permutation array of Bio::FastaFormat entries, with a small change in the fragment order
def self.mini_mutate (fasta)
	i = rand(fasta.length)
	j = rand(fasta.length)
	a = fasta[i]
	b = fasta[j]
	fasta[i] = b
	fasta[j] = a
	return fasta
end

# Input 0: A permutation array of Bio::FastaFormat entries (fragment arrangement)
# Input 1: Array of all the outputs from get_snp_data method
# Input 2: NOT IMPORTANT, WILL REMOVE EVENTUALLY "same" if we wish to compare the reformed ratio to the constant RATIO, or any other argument if we wish to compare to a re-run of comparable_ratio function in comparable_ratio.R
# Output: A correlation value that is the fitness of the Input 0 permutation
def self.fitness (fasta, snp_data, same) 
	id_n_lengths = fasta_id_n_lengths(fasta)
	fasta_ids = id_n_lengths[0]
	fasta_lengths = id_n_lengths[1]
	snps_per_frag = snps_per_fasta_frag(snp_data[2], fasta) #array of no. of snps per frag in same order as fasta
	pos_n_info = get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) #get snp positions for each frag in array of arrays
	pos = pos_n_info[0]
	info = pos_n_info[1] # rearranged vcfs_info filed to permutation order
	actual_pos = total_pos(pos, fasta_lengths)
	het_hom_snps = het_hom(actual_pos, info)
	het = het_hom_snps[0]
	hom = het_hom_snps[1]
	myr = RinRuby.new(echo=false)
	myr.assign "het_snps", het
	myr.assign "hom_snps", hom
	myr.eval "source('~/fragmented_genome_with_snps/lib/comparable_ratio.R')"
	if same == "diff"
		ratio = myr.pull "comparable_ratio(1)"
	else
		ratio = RATIO
	end
	myr.assign "ratio", ratio
	myr.eval "correlation <- get_corr(het_snps, hom_snps, ratio)"
	if same == "figure"
		myr.assign "dataset", ARGV[0]
		myr.eval "real_ratio <- get_real_ratio(het_snps, hom_snps, ratio)"
		myr.eval "plot_distribution(real_ratio, dataset)" # plot of the real_ratio
	end
	correlation = myr.pull "correlation"
	myr.quit
	return correlation
end

# Input 0: Array of Bio::FastaFormat entries
# Input 1: Integer of the desired population size
# Output: Population - Array of size "size", where each element is a shuffled permutation of the input 0 array of Bio::FastaFormat entries (random permutations)
def self.initial_population(fasta, size)
	population = []
	size.times do
		chromosome = fasta.shuffle
		population << chromosome
	end
	return population
end

# Input 0: Population - array of arrays where each sub array is a permutation of the fragments (Bio::FastaFormat entries)
# Input 1: Array of all the outputs from get_snp_data method
# Input 2: Integer of the desired number of permutations to be selected for the next generation
# Output 1: Array of fittest selection of Input 0 population: each sub array has two elements, the fitness and the permutation (which is itself an array of fragments)
# Output 2: Integer of leftover fragments, to be taken from the multiplied selected population
def self.select(pop, snp_data, num)
	fits = []
	pop.each do |sol| #solution
		fits << fitness(sol, snp_data, "same")
	end
	pop_fits = fits.zip(pop).sort
	length = pop_fits.length
	sliced = pop_fits.reverse.each_slice(num).to_a
	pop_fits = sliced[0].reverse
	if sliced[-1].length != sliced[0].length # if there is a remainder slice
		leftover = sliced[-1].length
	else 
		leftover = 0
	end
	#pop_fits = pop_fits.reverse.each_slice(length/2).to_a[0].reverse # best half gets saved
	puts "Selected #{pop_fits.size} of #{length} permutations"
	return pop_fits, leftover
end

# Input 0: Array of fittest selection of previous population: each sub array has two elements, the fitness and the permutation (which is itself an array of fragments)
# Input 1: Integer of the desired population size (NOT USED)
# Input 2: Integer of the desired number of mutant permutations in the new population (this number of mutate and mini_mutate methods)
# Input 3: Integer of the desired number of the best permutations from the previous population, to be included in the new one
# Input 4: Integer of the desired number of randomly shuffled permutations in a new population
# Input 5: Integer of leftover fragments, to be taken from the multiplied selected population
# Output: New population of mutants, recombinants etc - array of arrays where each sub array is a permutation of the fragments (Bio::FastaFormat entries)
def self.new_population(pop_fits, size, mut_num, save, ran, select_num, leftover) # mut_num = no. of mutants, save = number saved; from best, ran = no. of random permutations
	x = (size-leftover)/select_num
	pop_fits = pop_fits*x
	if leftover != 0
		pop_fits = [pop_fits, pop_fits[-leftover..-1]].flatten(1) #add leftover number of frags (best)
	end
	pop_save = pop_fits.reverse.each_slice(save).to_a[0] # saving best "save" of permutations
	pop = []
	pop_save.each do |i|
		pop << i[1] # adding the permutations only, not the fitness score
	end
	for i in pop_fits[(mut_num*2)+save+ran..-1]
		x = rand(size-1)
		pop << recombine(i[1], pop_fits[x][1])
	end
	mut_num.times do
		pop << mutate(pop_fits[rand(pop_fits.length)][1]) # mutating randomly selected permutations from pop_fits
	end
	mut_num.times do
		pop << mini_mutate(pop_fits[-1][1]) # mini_mutating the best permutations
	end
	ran.times do
		pop << pop_fits[0][1].shuffle
	end
	puts "New population size = #{pop.size}, with #{pop.size-(mut_num*2)-save-ran} recombinants, #{mut_num} mutants, #{mut_num} mini_mutants, the #{save} best from the previous generation and #{ran} random permutations."
	return pop
end

# Input 0: A permutation array of Bio::FastaFormat entries
# Input 1: VCF file
# Input 2: Integer of the number of times the fitness method should be called with the permutation, and to divide the sum by to get an average
# Output: Average fitness correlation value
def self.average_fitness (fasta, vcf_file, num)
	snp_data = get_snp_data(vcf_file)
	fits = []
	num.times do
		fits<<fitness(fasta, snp_data, "diff")
	end
	worst = fits.sort[0]
	average = fits.inject(:+)/num
	puts "Worst #{worst}"
	puts "Average #{average}"
	return average
end

# Input 0: FASTA file
# Input 1: VCF file
# Input 2: Integer of desired number of generations - the number of times a new population is created from an old one
# Input 3: Integer of desired size of each population (array of arrays where each sub array is a permutation of the fragments (Bio::FastaFormat entries))
# Input 4: Integer of the desired number of mutant permutations in each new population (this number of mutate and mini_mutate methods)
# Input 5: Integer of the desired number of the best permutations from each population, to be included in the next one
# Input 6: Integer of the desired number of randomly shuffled permutations in each new population
# Input 7: Correctly ordered array of Bio::FastaFormat entries
# Output 1: A saved .txt file of the fragment identifiers, of a permutation with a fitness that suggests it is the correct order
# Output 2: A saved figure of the algorithm's performance
# Output 3: A saved figure of the best permuation's homozygous/heterozygous SNP density ratio across the genome, assuming the fragment permutation is correct
def self.evolve(fasta_file, vcf_file, gen, pop_size, select_num, mut_num, save, ran, ordered_fasta, figures)
	gen_fits = [] # array of the best fitness in each generation
	ordered_ids = fasta_id_n_lengths(ordered_fasta)[0]
	snp_data = get_snp_data(vcf_file) #array of vcf frag ids, snp positions (fragments with snps), hash of each frag from vcf with no. snps, array of info field
	puts "Original order correlation = #{fitness(ordered_fasta, snp_data, "same")}"
	puts
	puts "Gen 0"
	fasta = fasta_array(fasta_file) #array of fasta format fragments
	pop = initial_population(fasta, pop_size)
	pop_fits_n_leftover = select(pop, snp_data, select_num)
	pop_fits = pop_fits_n_leftover[0]
	leftover = pop_fits_n_leftover[1]
	best_perm = pop_fits[-1][1]
	puts "Best correlation = #{pop_fits[-1][0]}"
	gen_fits << pop_fits[-1][0]
	best_perm_ids = fasta_id_n_lengths(best_perm)[0]
	puts
	y=1
	z=1
	gen.times do
		puts "Gen #{y}"
		prev_best_arr = pop_fits[-1][1]
		prev_best_fit = pop_fits[-1][0]
		pop = new_population(pop_fits, pop_size, mut_num, save, ran, select_num, leftover)
		pop_fits_n_leftover = select(pop, snp_data, select_num)
		pop_fits = pop_fits_n_leftover[0]
		leftover = pop_fits_n_leftover[1]
		best_perm = pop_fits[-1][1]
		puts "Best correlation = #{pop_fits[-1][0]}"
		gen_fits << pop_fits[-1][0]
		best_perm_ids = fasta_id_n_lengths(best_perm)[0]
		if pop_fits[-1][0] <= prev_best_fit
			puts "No fitness improvement" # If this is not called, this implies there has been some improvement
			z+=1
		else
			puts "FITNESS IMPROVEMENT!"
			z=1
		end
		puts
		if z >= 10
			puts "Algorithm quit for lack of fitness improvement, rearrangement score of #{rearrangement_score(ordered_ids, best_perm_ids)}"
			write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/reformed_ratio_frag_order.txt', best_perm_ids)
			fitness(best_perm, snp_data, "figure") # makes figure of ratio density distribution
		end
		if pop_fits[-1][0] >= 0.995 # If it looks like we have a winner, IN THE FINISHED ALGORITHM, THIS SHOULD BE...
			av = average_fitness(pop_fits[-1][1], vcf_file, 10)
			if av >= 0.999 && figures != "no figures"
				puts "Fitness #{av}: reform ratio has a rearrangement score of #{rearrangement_score(ordered_ids, best_perm_ids)}"
				write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/reformed_ratio_frag_order.txt', best_perm_ids)
				fitness(best_perm, snp_data, "figure") # makes figure of ratio density distribution
				z = 10
			end
		end
		if z >= 10
			then break
		end
		y+=1
		Signal.trap("PIPE", "EXIT")
	end
	if figures != "no figures"
		myr = RinRuby.new(echo=false)
		myr.assign "gen_fits", gen_fits
		myr.eval "source('~/fragmented_genome_with_snps/lib/comparable_ratio.R')"
		myr.assign "dataset", ARGV[0]
		myr.eval "plot_performance(gen_fits, dataset)"
		myr.quit
	end
end

end




