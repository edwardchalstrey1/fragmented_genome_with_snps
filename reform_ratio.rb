require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'

myr = RinRuby.new(echo=false)
myr.eval "source('~/fragmented_genome_with_snps/comparable_ratio.R')"
RATIO = myr.pull "comparable_ratio(1)"
myr.quit

def fasta_id_n_lengths (fasta) # array of fasta format frags
	ids = []
	lengths = []
	fasta.each do |i|
		ids << i.entry_id
		lengths << i.length
	end
	return ids, lengths
end
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
def get_positions (fasta, vcfs_chrom, vcfs_pos, snps_per_frag, vcfs_info)
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
def total_pos (pos, fasta_lengths)  # both args in same order as fasta = good. this all works!
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
	return totals
end
def het_hom (actual_pos, vcfs_info) #actual_pos in same order as fasta perm. now so is info
	het = []
	hom = []
	x = 0
	actual_pos.flatten.each do |snp|
		if vcfs_info.flatten[x] == {"AF"=>"1.0"} # homozygous SNPs have AF= 1.0, we can change this to a range for real data
			hom << snp
		elsif vcfs_info.flatten[x] == {"AF"=>"0.5"}
			het << snp
		end
		x+=1
	end
	return het, hom
end
def rearrangement_score (frags_original_order, rearranged)
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
def write_txt (filename, array)
	File.open(filename, "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end
# Genetic Algorithm
#------------------

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
			redo
		elsif x == mum.length # to compensate for datasets with a prime number of fragments:
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
		kid << child
	end
	return kid[0]
end
def mutate (fasta)
	x = 0
	until x > 2
		if prime?(fasta.length)
			x = division(fasta, "p")
		else
			x = division(fasta, "n")
		end
	end
	sliced = fasta.each_slice(x).to_a
	e = rand(sliced.length-1).to_i
	sliced[e] = sliced[e].shuffle
	return sliced.flatten
end
def mini_mutate (fasta)
	i = rand(fasta.length)
	j = rand(fasta.length)
	a = fasta[i]
	b = fasta[j]
	fasta[i] = b
	fasta[j] = a
	return fasta
end
def fitness (fasta, snp_data, same) # same? - use the same constant RATIO of snp ratio to qq plot against the frags, or not
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
	myr.eval "source('~/fragmented_genome_with_snps/comparable_ratio.R')"
	if same == "same"
		ratio = RATIO
	else
		ratio = myr.pull "comparable_ratio(1)"
	end
	myr.assign "ratio", ratio
	coeff = myr.pull "qq_real_expect(het_snps, hom_snps, ratio)"
	myr.quit
	return coeff
end
def initial_population(fasta, size)
	population = []
	size.times do
		chromosome = fasta.shuffle
		population << chromosome
	end
	return population
end
def select(pop, snp_data)
	fits = []
	pop.each do |sol| #solution
		fits << fitness(sol, snp_data, "same")
	end
	pop_fits = fits.zip(pop).sort
	length = pop_fits.length
	pop_fits = pop_fits.reverse.each_slice(length/2).to_a[0].reverse # best half gets saved
	puts "Selected #{pop_fits.size} of #{length} permutations"
	return pop_fits
end
def new_population(population, size, mut_num, save, ran) # mut_num = no. of mutants, save = number saved; from best, ran = no. of random permutations
	population = [population, population].flatten(1)
	pop = []
	population[-save,save].each do |i|
		pop << i[1]
	end
	x = rand(size-1)
	for i in population[(mut_num*2)+save+ran..-1]
		pop << recombine(i[1], population[x][1])
	end
	mut_num.times do
		pop << mutate(population[rand(population.length)][1])
	end
	mut_num.times do
		pop << mini_mutate(population[-1][1])
	end
	ran.times do
		pop << population[0][1].shuffle
	end
	return pop
end
def average_fitness (fasta, vcf_file, num)
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
def evolve(fasta_file, vcf_file, gen, pop_size, mut_num, save, ran, ordered_fasta)
	ordered_ids = fasta_id_n_lengths(ordered_fasta)[0]
	snp_data = get_snp_data(vcf_file) #array of vcf frag ids, snp positions (fragments with snps), hash of each frag from vcf with no. snps, array of info field
	puts "Original order correlation = #{fitness(ordered_fasta, snp_data, "same")}"
	#puts "Original order score = #{rearrangement_score(ordered_ids, ordered_ids)}"
	#worst_score = rearrangement_score(ordered_ids, ordered_ids.reverse)
	puts
	puts "Gen 0"
	fasta = fasta_array(fasta_file) #array of fasta format fragments
	pop = initial_population(fasta, pop_size)
	pop_fits = select(pop, snp_data)
	puts "Best correlation = #{pop_fits[-1][0]}"
	best_perm_ids = fasta_id_n_lengths(pop_fits[-1][1])[0]
	#puts "Score of best correlation = #{rearrangement_score(ordered_ids, best_perm_ids)}  Worst score = #{worst_score}" # THE REARRANGEMENT SCORE METHOD IS FOR TESTING THE ALGORITHM ONLY, NOT PART OF IT
	puts
	y=1
	z=1
	gen.times do
		puts "Gen #{y}"
		prev_best_arr = pop_fits[-1][1]
		prev_best_fit = pop_fits[-1][0]
		pop = new_population(pop_fits, pop_size, mut_num, save, ran)
		pop_fits = select(pop, snp_data)
		puts "Best correlation = #{pop_fits[-1][0]}"
		best_perm_ids = fasta_id_n_lengths(pop_fits[-1][1])[0]
		#score = rearrangement_score(ordered_ids, best_perm_ids)
		#puts "Score of best correlation = #{score}  Worst score = #{worst_score}" # THE REARRANGEMENT SCORE METHOD IS FOR TESTING THE ALGORITHM ONLY, NOT PART OF IT
		if pop_fits[-1][0] <= prev_best_fit
			puts "No fitness improvement" # If this is not called, this implies there has been some improvement
			z+=1
		else
			puts "FITNESS IMPROVEMENT!"
			z=1
		end
		#if score < rearrangement_score(ordered_ids, (fasta_id_n_lengths(prev_best_arr)[0]))
		#	puts "SCORE IMPROVEMENT!"
		#else
		#	puts "No score improvement"
		#end
		puts
		if z >= 10
			puts best_perm_ids
		end
		if z >= 10
			then break
		end
		if pop_fits[-1][0] >= 0.995 # If it looks like we have a winner, IN THE FINISHED ALGORITHM, THIS SHOULD BE...
			av = average_fitness(pop_fits[-1][1], vcf_file, 10)
			if av >= 0.999
				puts "Perfect reform ratio has a rearrangement score of #{rearrangement_score(ordered_ids, best_perm_ids)}"
				write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/reformed_ratio_frag_order.txt', best_perm_ids)		
			end
		end
		y+=1
		Signal.trap("PIPE", "EXIT")
	end
end

vcf = 'arabidopsis_datasets/'+ARGV[0].to_s+'/snps.vcf'
fasta = 'arabidopsis_datasets/'+ARGV[0].to_s+'/frags_shuffled.fasta'

ordered_fasta = fasta_array('arabidopsis_datasets/'+ARGV[0].to_s+'/frags.fasta')
#average_fitness(ordered_fasta, vcf, 10) # test to see how well correct arrangement performs...
#average_fitness(fasta_array(fasta), vcf, 10) # ... vs random arrangement

evolve(fasta, vcf, 200, 100, 10, 10, 5, ordered_fasta) # gen, pop, mut*2, save, ran ### ordered_ids is temporary








def how_we_doin (ordered_fasta, vcf)
	ordered_ids = fasta_id_n_lengths(ordered_fasta)[0]
	snp_data = get_snp_data(vcf)
	puts "Ordered #{fitness(ordered_fasta, snp_data, "same")}"
	puts rearrangement_score(ordered_ids, ordered_ids)
	mutant = mutate(ordered_fasta)
	puts "Mutant #{fitness(mutant, snp_data, "same")}"
	puts rearrangement_score(ordered_ids, fasta_id_n_lengths(mutant)[0])
	recombination = recombine(ordered_fasta, ordered_fasta.shuffle)
	puts "recombination #{fitness(recombination, snp_data, "same")}"
	puts rearrangement_score(ordered_ids, fasta_id_n_lengths(recombination)[0])
	#fits = []
	#100.times {|x| fits << fitness(ordered_fasta.shuffle, snp_data, "same")}
	#fits = fits.inject(:+)/100
	shuffled = ordered_fasta.shuffle
	puts "Shuffled #{fitness(shuffled, snp_data, "same")}"
	puts rearrangement_score(ordered_ids, fasta_id_n_lengths(shuffled)[0])
	#puts "Shuffled #{fits}"
	mini = mini_mutate(ordered_fasta)
	puts "mini mutate #{fitness(mini, snp_data, "same")}"
	puts rearrangement_score(ordered_ids, fasta_id_n_lengths(mini)[0])
end
#how_we_doin(ordered_fasta, vcf)


