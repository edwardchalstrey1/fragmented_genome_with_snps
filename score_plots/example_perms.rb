class ExamplePerms

	require '~/fragmented_genome_with_snps/lib/reform_ratio.rb'
	require '~/fragmented_genome_with_snps/lib/GATOC.rb'

	# Input 0: Array of fasta format fragments
	# Input 1: Permutation: txt file saved with frag ids
	# Output: Array of the fitness for the permutation with the fatsa frags re-ordered by the permutation
	def self.fasta_p(fasta, perm)
		perm_ids = []
		IO.foreach(perm) { |line| perm_ids << line.gsub(/\n/,'') }
		perm_ids = perm_ids[1..-1]
		fitness = perm_ids[0]
		ids = ReformRatio::fasta_id_n_lengths(fasta)[0]
		fasta_perm = []
		perm_ids.each do |id|
			fasta.each do |frag|
				if id == frag.entry_id
					fasta_perm << frag
				end
			end
		end
		return [fitness, fasta_perm].flatten
	end

	# Input 0: Array of fasta format fragments
	# Input 1: Integer of the number of permutations desired for each population
	# Output: Array of populations, each population contains permutations of a certain type. Permutations are arrays of fasta frags, with fitness
	def self.get_perms(fasta, pop_size, snp_data)
		mut_pop, shuf_pop = [], []
		pop_size.times do
			mut_perm = GATOC::mutate(fasta)
			fitness = GATOC::fitness(mut_perm, snp_data, 'same')
			mut_perm_id = [fitness, ReformRatio::fasta_id_n_lengths(mut_perm)[0]].flatten
			mut_pop << mut_perm_id

			shuf_perm = fasta.shuffle
			fitness = GATOC::fitness(shuf_perm, snp_data, 'same')
			shuf_perm_id = [fitness, ReformRatio::fasta_id_n_lengths(shuf_perm)[0]].flatten
			shuf_pop << shuf_perm_id
		end
		return mut_pop, shuf_pop
	end
end


# best_g10 = fasta_p(fasta, "~/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen10/best_permutation.txt")