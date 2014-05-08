#encoding: utf-8
class GATOC # Genetic Algorithm To Order Contigs
	require_relative 'snp_dist'
	require_relative 'write_it'
	require_relative 'reform_ratio'
	require '~/pmeth/lib/pmeth'

	# Input 0: A permutation array of Bio::FastaFormat entries (fragment arrangement)
	# Input 1: Array of all the outputs from get_snp_data method
	# Input 2: Integer of the generation that the genetic algorithm is on, if gen is not an integer, no figure is plotted
	# Input 3: Example ratio to compare against in Q-Q plot
	# Input 4: Location to save plots e.g. "~/fragmented_genome_with_snps/arabidopsis_datasets"
	# Input 5: Dataset algorithm running on
	# Input 6: Name of this run of the algorithm
	# Input 7: Length of divisions of the genome to calculate the SNP frequency of
	# Input 8: Length of the genome
	# Output 0: A correlation value that is the fitness of the Input 0 permutation
	# Output 1: List of homozygous SNPs for the permutation
	# Output 2: Heterozygous list
	# Output 3: List aof values representing the ratio distribution
	def self.fitness(fasta, snp_data, gen, comparable_ratio, location, dataset, run, div, genome_length)
		snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) #array of no. of snps per frag in same order as fasta
		pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) #get snp positions for each frag in array of arrays
		actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
		het_snps, hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])
		fratio_breaks_perm = SNPdist::fratio(hom_snps, het_snps, div, genome_length) # frequency ratio array
		perm_ratio = SNPdist::hyp_snps(fratio_breaks_perm, div, genome_length) # hypothetical snp positions array
		correlation = SNPdist::qq_cor(comparable_ratio, perm_ratio)
		if Integer === gen # create figure of the hypothetical snp distribution over the genome, for each generation
			SNPdist::plot_hyp(perm_ratio, location, "#{dataset}/#{run}", gen, genome_length)
		end
		return correlation, hom_snps, het_snps, perm_ratio
	end

	# Input 0: Array of Bio::FastaFormat entries (or any array)
	# Input 1: Integer of the desired population size
	# Output: Population - Array of size "size", where each element is a shuffled permutation of the input 0 array of Bio::FastaFormat entries (random permutations)
	def self.initial_population(fasta, size)
		population = []
		size.times do
			population << [fasta.shuffle, 'random']
		end
		return population
	end

	# Input 0: Population - array of arrays where each sub array is a permutation of the fragments (Bio::FastaFormat entries)
	# Input 1: Array of all the outputs from get_snp_data method
	# Input 2: Integer of the desired number of permutations to be selected for the next generation
	# Input 3: The example homozygous/heterozygous SNP ratio to compare against in fitness
	# Input 4: Length of divisions of the genome to calculate the SNP frequency of
	# Input 5: Length of the genome
	# Output 0: Array of fittest selection of Input 0 population: each sub array has two elements, the fitness and the permutation (which is itself an array of fragments)
	# Output 1: Integer of leftover permutations, to be taken from the multiplied selected population
	# Output 2: Pre-selected version of output 0
	# Output 3: Array of strings, each string represents the method used to create a permutation followed by fitness score, and the array is ordered the same as output 2
	def self.select(pop, snp_data, num, ratio, div, genome_length)
		puts "Pop is unique: #{pop.uniq.length == pop.length}"
		fits = {}
		pop.each do |fasta_array, type|
			fitn = fitness(fasta_array, snp_data, 'same', ratio, 'location not needed', 'dataset not needed', 'run not needed', div, genome_length)[0]
			fits[fasta_array] = [fitn, type] # maybe some have exact same fitness, perhaps we can make fitness the value, then sort by value
		end
		if fits.size < pop.size # to compensate for duplicates, we add extra shuffled permutations #TODO change to swap mutants
			diff = pop.size - fits.size
			x = 0
			diff.times do
				extra_rand = pop[0][0].shuffle
				fitn = fitness(extra_rand, snp_data, 'same', ratio, 'location not needed', 'dataset not needed', 'run not needed', div, genome_length)[0]
				fits[extra_rand] = [fitn, 'extra_rand']
				x+=1
			end
			puts "#{x} extra random permutations added, due to multiples of the same permutation in the population"
		end
		fits = fits.sort_by {|k,v| v[0]} # sorting by fitness score
		types = []
		fits.each {|k,v| types << v.reverse} # adding the types with fitness to a new array
		x = 0
		fits.each {|k,v| fits[x][1] = v[0]; x+=1} # getting rid of the types, so v is now just fitness score
		pop_fits = []
		fits.each {|i| pop_fits << i.reverse} # swapping the permutation/fitness score around
		initial_pf = pop_fits # the input permutations ordered by fitness, not yet selcted
		sliced = pop_fits.reverse.each_slice(num).to_a # sliced the population ordered by fitness into chunks of size num, choosing the chunk with the highest fitness scores (reversing and choosing chunk 0)
		pop_fits = sliced[0].reverse # creating the selected population, and reversing them to ascending fitness order
		if sliced[-1].length != sliced[0].length # if there is a remainder slice
			leftover = sliced[-1].length
		else 
			leftover = 0
		end
		return pop_fits, leftover, initial_pf, types.flatten # flatten types to get type then fitness alternate
	end

	# Input 0: Array of fittest selection of previous population: each sub array has two elements, the fitness and the permutation (which is itself an array of fragments)
	# Input 1: Integer of the desired population size
	# Input 2: Integer of the desired number of chunk mutant permutations in the new population
	# Input 3: Integer of the desired number of swap mutant permutations in the new population
	# Input 4: Integer of the desired number of the best permutations from the previous population, to be included in the new one
	# Input 5: Integer of the desired number of randomly shuffled permutations in a new population
	# Input 6: Integer of the number of permutations selected by the select method
	# Input 7: Integer of leftover permutations, to be taken from the multiplied selected population
	# Output: New population - array of arrays where each sub array is a permutation of the fragments (Bio::FastaFormat entries)
	def self.new_population(pop_fits, size, c_mut, s_mut, save, ran, select_num, leftover) # mut_num = no. of mutants, save = number saved; from best, ran = no. of random permutations
		x = (size-leftover)/select_num
		pop_fits = pop_fits * x
		if leftover != 0
			pop_fits = [pop_fits, pop_fits[-leftover..-1]].flatten(1) # add leftover number of frags (best)
			puts "#{leftover} leftover frags added"
		end
		pop_save = pop_fits.reverse.each_slice(save).to_a[0] # saving best "save" of permutations # TODO are types wrong?
		pop = []
		pop_save.each{|i| pop << [i[1], 'saved']} # adding the permutations only, not the fitness score
		c_mut.times{pop << [PMeth.chunk_mutate(pop_fits[rand(pop_fits.length)][1].dup), 'chunk_mutant']} # chunk_mutating randomly selected permutations from pop_fits
		s_mut.times{pop << [PMeth.swap_mutate(pop_fits[-1][1].dup), 'swap_mutant']} # swap_mutating the best permutations
		ran.times{pop << [pop_fits[0][1].shuffle, 'random']}
		new_pop_msg = "Population size = #{pop.size}, with #{c_mut} chunk mutants, #{s_mut} swap mutants, the #{save} best from the previous generation and #{ran} random permutations."
		return pop, new_pop_msg
	end

	# Input 0: Fittest selection array: see output 0 for select method
	# Input 1: Location to save files to e.g. "fragmented_genome_with_snps/arabidopsis_datasets"
	# Input 2: Dataset algorithm running on
	# Input 3: Name of this run of the algorithm
	# Input 4: Generation of the genetic algorithm
	def self.save_perms(pop_fits, location, dataset, run, gen)
		Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}/Gen#{gen}"))
		Dir.chdir(File.join(Dir.home, "#{location}")) do
			x = 1
			pop_fits.each do |fitness, perm|
				ids = ReformRatio::fasta_id_n_lengths(perm)[0]
				WriteIt::write_txt("#{dataset}/#{run}/Gen#{gen}/permutation#{x}", [fitness, ids].flatten)
				x+=1
			end
			if gen != 0
				ids = ReformRatio::fasta_id_n_lengths(pop_fits[-1][1])[0]
				WriteIt::write_txt("#{dataset}/#{run}/Gen#{gen}/best_permutation", [pop_fits[-1][0], ids].flatten) # fitness and ids
			end
		end
	end

	# Input 0: FASTA file
	# Input 1: VCF file
	# Input 2: parameters:
	# 	gen: Integer of desired number of generations - the number of times a new population is created from an old one
	# 	pop_size: Integer of desired size of each population (array of arrays where each sub array is a permutation of the fragments (Bio::FastaFormat entries))
	# 	c_mut: Integer of the desired number of chunk mutant permutations in each new population
	#   s_mut: Integer of the desired number of swap mutant permutations in each new population
	# 	save: Integer of the desired number of the best permutations from each population, to be included in the next one
	# 	ran: Integer of the desired number of randomly shuffled permutations in each new population
	#   loc: Location to save output files to
	#   comparable_ratio: An example distribution to compare the ratio of homozygous/heterozygous SNP density to
	#   dataset: The sub folder containing fasta and vcf files
	#   run: The name you'd like to assign this run of the algorithm
	#   div: Length of divisions of the genome to calculate the SNP frequency of, in fitness method
	#   genome_length: Length of the genome the algorithm is being run on
	#   end: When the fitness score of a permutation generated by the algorithm reaches this value, the algorithm will stop
	#   gen_end: The number of generations the algorithm can go without creating a permutation with increased fitness, before quitting
	# Output 1: A saved .txt file of the fragment identifiers, of a permutation with a fitness that suggests it is the correct order
	# Output 2: A saved figure of the algorithm's performance
	# Output 3: A saved figure of the best permuation's homozygous/heterozygous SNP density ratio across the genome, assuming the fragment permutation is correct
	def self.evolve(fasta_file, vcf_file, parameters)
		opts = {
			:gen => 10000000000,
			:pop_size => 100,
			:select_num => 50,
			:c_mut => 10,
			:s_mut => 10,
			:save => 5,
			:ran => 5,
			:loc => '~/fragmented_genome_with_snps/arabidopsis_datasets',
			:comparable_ratio => 'nothing', # TODO add in default ratio
			:dataset => ARGV[0],
			:run => ARGV[1],
			:div => 100.0,
			:genome_length => 2000.0,
			:end => 0.999,
			:gen_end => 10
			}.merge!(parameters)
		# Dir.mkdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}")) # make the directory to put data files into

		gen_fits = [] # array of the best fitness in each generation
		snp_data = ReformRatio::get_snp_data(vcf_file) #array of vcf frag ids, snp positions (fragments with snps), hash of each frag from vcf with no. snps, array of info field
		fasta = ReformRatio::fasta_array(fasta_file) #array of fasta format fragments
		pop = initial_population(fasta, opts[:pop_size])

		pop_fits, leftover, initial_pf, types = select(pop, snp_data, opts[:select_num], opts[:comparable_ratio], opts[:div], opts[:genome_length])
		fit, hm, ht, hyp = fitness(pop_fits[-1][1], snp_data, 0, opts[:comparable_ratio], opts[:loc], opts[:dataset], opts[:run], opts[:div], opts[:genome_length]) # makes figure of ratio density distribution for the best permutation in each generation
		Dir.mkdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}/Gen0_lists"))
		Dir.chdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}/Gen0_lists")) do
			WriteIt::write_txt("gen_0_hm", hm)
			WriteIt::write_txt("gen_0_ht", ht)
			WriteIt::write_txt("gen_0_hyp", hyp)
			WriteIt::write_txt("gen_0_types", types)
		end
		puts "Gen0 \n Best correlation = #{pop_fits[-1][0]}\n \n"
		gen_fits << pop_fits[-1][0]
		save_perms(initial_pf, opts[:loc], opts[:dataset], opts[:run], 0)

		gen, z, messages = 1, 1, [] # can check messages for errors
		opts[:gen].times do

			prev_best_fit = pop_fits[-1][0]
			pop, msg = new_population(pop_fits, opts[:pop_size], opts[:c_mut], opts[:s_mut], opts[:save], opts[:ran], opts[:select_num], leftover)
			messages << msg
			pop_fits, leftover, initial_pf, types = select(pop, snp_data, opts[:select_num], opts[:comparable_ratio], opts[:div], opts[:genome_length])
			gen_fits << pop_fits[-1][0]
			save_perms(initial_pf, opts[:loc], opts[:dataset], opts[:run], gen)

			puts "Gen#{gen}\n Best correlation = #{pop_fits[-1][0]}"
			if pop_fits[-1][0] <= prev_best_fit
				puts "No fitness improvement\n \n" # If this is not called, this implies there has been some improvement
				z+=1
			else
				puts "FITNESS IMPROVEMENT!\n \n"
				z = 1
			end
				
			fit, hm, ht, hyp = fitness(pop_fits[-1][1], snp_data, gen, opts[:comparable_ratio], opts[:loc], opts[:dataset], opts[:run], opts[:div], opts[:genome_length]) # makes figure of ratio densitgen distribution for the best permutation in each generation
			Dir.mkdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}/Gen#{gen}_lists"))
			Dir.chdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}/Gen#{gen}_lists")) do
				WriteIt::write_txt("gen_#{gen}_hm", hm)
				WriteIt::write_txt("gen_#{gen}_ht", ht)
				WriteIt::write_txt("gen_#{gen}_hyp", hyp)
				WriteIt::write_txt("gen_#{gen}_types", types)
			end

			if z >= opts[:gen_end] || gen >= opts[:gen] || pop_fits[-1][0] >= opts[:end]
				then break
			end
			gen+=1
			Signal.trap("PIPE", "EXIT")
		end
		Dir.chdir(File.join(Dir.home, opts[:loc])) do
			WriteIt::write_txt("#{opts[:dataset]}/#{opts[:run]}/reformed_ratio_frag_order", 
				[opts, "generation #{gen}/#{opts[:gen]}: #{pop_fits[-1][0]}", pop_fits[-1][0],'###', 
				ReformRatio::fasta_id_n_lengths(pop_fits[-1][1])[0]].flatten)
			WriteIt::write_txt("#{opts[:dataset]}/#{opts[:run]}/messages", messages)
		end
	end
end