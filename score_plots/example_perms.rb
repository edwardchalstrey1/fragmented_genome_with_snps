class ExamplePerms

	require '~/fragmented_genome_with_snps/lib/reform_ratio.rb'
	require '~/fragmented_genome_with_snps/lib/GATOC.rb'

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
		return fasta_perm, fitness
	end

	def self.get_perms(fasta, pop_size)
		mut_pop, shuf_pop = [], []
		pop_size.times do
			mut_pop << ReformRatio::fasta_id_n_lengths(GATOC::mutate(fasta))[0]
			shuf_pop << ReformRatio::fasta_id_n_lengths(fasta.shuffle)[0]
		end
		return mut_pop, shuf_pop
	end
end


# best_g10 = fasta_p(fasta, "~/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen10/best_permutation.txt")