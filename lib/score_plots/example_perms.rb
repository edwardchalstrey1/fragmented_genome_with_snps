#encoding: utf-8
class ExamplePerms

	require_relative '../reform_ratio.rb'
	require_relative '../GATOC.rb'
	require 'pmeth'

	# Input 0: Array of fasta format fragments
	# Input 1: Permutation: txt file saved with frag ids
	# Output: Array of the fitness for the permutation with the fatsa frags re-ordered by the permutation
	def self.fasta_p(fasta, perm)
		perm_ids = []
		IO.foreach(perm) { |line| perm_ids << line.gsub(/\n/,'') }
		perm_ids = perm_ids[1..-1]
		fitness = perm_ids[0].to_f
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
	# Input 2: Output from ReformRatio::get_snp_data on VCF file
	# Input 3: Ratio vector (array) of distribution for correctly ordered contigs
	# Input 4: Number of divisions at which to calculate SNP density in permutation
	# Input 5: Genome length
	# Output: Array of populations, each population contains permutations of a certain type. Permutations are arrays of fasta frags, with fitness
	def self.get_perms(fasta, pop_size, snp_data, comparable_ratio, div, genome_length)
		chunk, swap, shuf = [], [], []
		pop_size.times do
			chunk_perm = PMeth.chunk_mutate(fasta.dup)
			fitness = GATOC::fitness(chunk_perm, snp_data, 'same', comparable_ratio, 'location', 'dataset', 'run', div, genome_length)[0]
			chunk_perm_id = [fitness, ReformRatio::fasta_id_n_lengths(chunk_perm)[0]].flatten
			chunk << chunk_perm_id

			swap_perm = PMeth.swap_mutate(fasta.dup)
			fitness = GATOC::fitness(swap_perm, snp_data, 'same', comparable_ratio, 'location', 'dataset', 'run', div, genome_length)[0]
			swap_perm_id = [fitness, ReformRatio::fasta_id_n_lengths(swap_perm)[0]].flatten
			swap << swap_perm_id

			shuf_perm = fasta.shuffle
			fitness = GATOC::fitness(shuf_perm, snp_data, 'same', comparable_ratio, 'location', 'dataset', 'run', div, genome_length)[0]
			shuf_perm_id = [fitness, ReformRatio::fasta_id_n_lengths(shuf_perm)[0]].flatten
			shuf << shuf_perm_id
		end
		return chunk, swap, shuf
	end
end

