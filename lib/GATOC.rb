#encoding: utf-8
class GATOC # Genetic Algorithm To Order Contigs
	require_relative 'fitness_score'
	require_relative 'write_it'
	require_relative 'reform_ratio'
	require 'pmeth'
	require_relative 'quit_if'

	# Input 0: A permutation array of Bio::FastaFormat entries (fragment arrangement)
	# Input 1: Array of all the outputs from get_snp_data method
	# Input 2: Example ratio to compare against in Q-Q plot
	# Input 3: Length of divisions of the genome to calculate the SNP frequency of
	# Input 4: Length of the genome
	# Output 0: A correlation value that is the fitness of the Input 0 permutation
	# Output 1: List of homozygous SNPs for the permutation
	# Output 2: Heterozygous list
	# Output 3: List of values representing the ratio distribution
	def self.fitness(fasta, snp_data, comparable_ratio, div, genome_length)
		het_snps, hom_snps = ReformRatio.perm_pos(fasta, snp_data)
		perm_ratio = FitnessScore.ratio(hom_snps, het_snps, div, genome_length)
		# correlation = FitnessScore::score(comparable_ratio, perm_ratio)
		correlation = FitnessScore.distance_score(hom_snps)
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
		fits = {}
		pop.each do |fasta_array, type|
			fitn = fitness(fasta_array, snp_data, ratio, div, genome_length)[0]
			fits[fasta_array] = [fitn, type] # maybe some have exact same fitness, perhaps we can make fitness the value, then sort by value
		end
		if fits.size < pop.size # to compensate for duplicates, we add extra swap mutants
			diff = pop.size - fits.size
			x = 0
			diff.times do
				extra_swap = PMeth.swap_mutate(pop[rand(pop[0].length)][0].dup)
				fitn = fitness(extra_swap, snp_data, ratio, div, genome_length)[0]
				fits[extra_swap] = [fitn, 'extra_swap']
				x+=1
			end
			puts "#{x} extra swap mutants added, due to multiples of the same permutation in the population"
		end
		fits = fits.sort_by {|k,v| v[0]} # sorting by fitness score
		types = []
		fits.each {|k,v| types << v.reverse} # adding the types with fitness to a new array
		x = 0
		fits.each {|k,v| fits[x][1] = v[0]; x+=1} # getting rid of the types, so v is now just fitness score
		pop_fits = []
		fits.each {|i| pop_fits << i.reverse} # swapping the permutation/fitness score around
		initial_pf = pop_fits.reverse # the input permutations ordered by fitness, not yet selcted ### (lowest and best score is last)
		sliced = pop_fits.each_slice(num).to_a # sliced the population ordered by fitness into chunks of size num, choosing the chunk with the highest fitness scores (reversing and choosing chunk 0)
		pop_fits = sliced[0].reverse # creating the selected population, and reversing them to ascending fitness order ### (lowest and best score is last)
		if sliced[-1].length != sliced[0].length # if there is a remainder slice
			leftover = sliced[-1].length
		else 
			leftover = 0
		end
		return pop_fits, leftover, initial_pf, types.reverse
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
		pop_save = pop_fits.reverse.each_slice(save).to_a[0] # saving best "save" of permutations
		pop = []
		pop_save.each{|i| pop << [i[1], 'saved']} # adding the permutations only, not the fitness score
		c_mut.times{pop << [PMeth.chunk_mutate(pop_fits[rand(pop_fits.length)][1].dup), 'chunk_mutant']} # chunk_mutating randomly selected permutations from pop_fits
		s_mut.times{pop << [PMeth.swap_mutate(pop_fits[-1][1].dup), 'swap_mutant']} # swap_mutating the best permutations
		ran.times{pop << [pop_fits[0][1].shuffle, 'random']}
		return pop
	end

	# Input 0: See output 2 of select
	# Input 1: Location to save files to
	# Input 2: Dataset algorithm running on
	# Input 3: Name of this run of the algorithm
	# Input 4: Generation of the genetic algorithm
	# Input 5: Output 3 of select
	# Output: txt files with data interpretable by text-table gem AND txt files for each permutation
	def self.save_perms(pop_fits, location, dataset, run, gen, types)
		Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}/Gen#{gen}"))
		Dir.chdir(File.join(Dir.home, "#{location}")) do
			table_data = [['Permutation', 'Fitness Score', 'Type', 'FASTA ids']]
			x = 1
			pop_fits.each do |fitness, permutation|
				table_data << ["permutation#{x}", fitness, types[x-1][0], ReformRatio::fasta_id_n_lengths(permutation)[0].join(", ")]
				x+=1
			end
			WriteIt::write_txt("#{dataset}/#{run}/Gen#{gen}/table_data", table_data)

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
	#   start_pop: Population array of permutations (arrays of FASTA format fragments)
	#   start_gen: Generation that the algorithm ist starting from
	#   auc: quit algorithm if area under curve of fitness improvement is < auc% better for this auc_gen generations than previous auc_gen
	#   auc_gen: see above
	#   restart_zero: should be set to anything other than nil, if the algorithm needs to be restarted from a generation zero population
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
			:comparable_ratio => nil,
			:dataset => ARGV[0],
			:run => ARGV[1],
			:div => 100.0,
			:genome_length => 2000.0,
			:start_pop => nil,
			:start_gen => 0,
			:auc => 5,
			:auc_gen => 5,
			:restart_zero => nil
			}.merge!(parameters)

		snp_data = ReformRatio::get_snp_data(vcf_file) # array of vcf frag ids, snp positions (fragments with snps), hash of each frag from vcf with no. snps, array of info field
		fasta = ReformRatio::fasta_array(fasta_file) # array of fasta format fragments
		prev_best_fit = nil
		if opts[:start_pop] == nil
			pop = initial_population(fasta, opts[:pop_size])
		else
			pop = opts[:start_pop]
		end

		gen, last_best, auc = opts[:start_gen], [], nil
		opts[:gen].times do

			if opts[:start_gen] == 0 && opts[:restart_zero] != nil && gen == 0
				gen+=1
			end

			pop_fits, leftover, initial_pf, types = select(pop, snp_data, opts[:select_num], opts[:comparable_ratio], opts[:div], opts[:genome_length])

			unless opts[:start_pop] != nil && gen == opts[:start_gen] # if using a starting population, we don't want to overwite files for that generation
				save_perms(initial_pf, opts[:loc], opts[:dataset], opts[:run], gen, types)
			end

			puts "Gen#{gen}\n Best correlation = #{pop_fits[-1][0]}"
			if prev_best_fit != nil
				if pop_fits[-1][0] >= prev_best_fit
					puts "No fitness improvement\n \n"
				else
					puts "FITNESS IMPROVEMENT!\n \n"
				end
			end
				
			fit, hm, ht, ratios = fitness(pop_fits[-1][1], snp_data, opts[:comparable_ratio], opts[:div], opts[:genome_length])

			unless opts[:start_pop] != nil && gen == opts[:start_gen] # if using a starting population, we don't want to overwite files for that generation
				Dir.mkdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}/Gen#{gen}_lists"))
				Dir.chdir(File.join(Dir.home, "#{opts[:loc]}/#{opts[:dataset]}/#{opts[:run]}/Gen#{gen}_lists")) do
					WriteIt::write_txt("gen_#{gen}_hm", hm)
					WriteIt::write_txt("gen_#{gen}_ht", ht)
					WriteIt::write_txt("gen_#{gen}_ratios", ratios)
				end
			end

			prev_best_fit = pop_fits[-1][0]
			pop = new_population(pop_fits, opts[:pop_size], opts[:c_mut], opts[:s_mut], opts[:save], opts[:ran], opts[:select_num], leftover)
			gen+=1

			last_best << prev_best_fit
			if last_best.length == opts[:auc_gen]
				last_auc = auc
				auc = QuitIf::quit(last_best)
				if last_auc != nil && (last_auc - auc) <= (last_auc/(100.0/opts[:auc].to_f))
					puts 'auc break'
					gen = opts[:gen]
				end
				last_best = []
			end

			# if pop_fits[-1][0] == 1.0
			# 	puts 'correct permutation achieved'
			# 	gen = opts[:gen]
			# end
			
			if gen >= opts[:gen]
				then break
			end
			Signal.trap("PIPE", "EXIT")
		end
	end
end