#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'
require_relative 'lib/write_it'
require_relative 'lib/example_perms'

dataset = ARGV[0]
run = ARGV[1]
gen =  ARGV[2].to_i # Number of generations
pop_size = ARGV[3].to_i # Population size
select_num = ARGV[4].to_i # Number of permutations to select from each generation
c_mut = ARGV[5].to_i # Number of chunk mutants (see GATOC and PMeth gem) in each new population
s_mut = ARGV[6].to_i # Number of swap mutants (see GATOC and PMeth gem) in each new population
save = ARGV[7].to_i # Number of permutations to save from each generation
ran = ARGV[8].to_i # Number of random permutations in each generation
restart = ARGV[9] # Tells the algorithm to continue from the most recent generation if it has stopped

## Files ##
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags.fasta"
location = 'fragmented_genome_with_snps/arabidopsis_datasets'

### RESTART ###
### Whilst testing the algorithm with my model genome, the code below works out what generation to restart from, if the algorithm stops unexpectedly
	pop, restart_gen, restart_zero = [], [], []
	fasta = ReformRatio::fasta_array(fasta_file)
	if restart == nil
		Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) # make the directory to put permutation files into WARNING: this NEEDS to be in the algorithm run script
		### Whilst testing the algorithm with my model genome, we need to know the correct permutation's SNP distributions
			snp_data = ReformRatio.get_snp_data(vcf_file)
			ht, hm = ReformRatio.perm_pos(fasta, snp_data)
			WriteIt::write_txt("arabidopsis_datasets/#{dataset}/hm_snps", hm)
			WriteIt::write_txt("arabidopsis_datasets/#{dataset}/ht_snps", ht)
		###
		restart_gen << 0
		pop = nil
	else ## For restarts ##
		Dir.chdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) do # Selecting the generation to restart from
			dirs = Dir.glob('*').select {|f| File.directory? f}
			dirs.delete("Gencorrect_lists")
			if dirs.length < 2
				restart_gen << 0
			elsif dirs.length == 2
				restart_zero <<'restart'
				restart_gen << 0
			else
				restart_gen << (dirs.length/2)-1
			end
		end
		if restart_gen[0] != 0
			id_pop = ExamplePerms::get_perms(1, restart_gen[0], 0, dataset, run).flatten(1) # should be population array of permutations (arrays of ids)
			id_pop.each do |ids|
				pop << [ExamplePerms::fasta_p_id(fasta, ids), 'type'] # move the 2 methods used here to a class that is not otherwise redundant
			end
		else
			pop = nil
		end
	end
###

### Run the algorithm ### (:start_pop should be set to nil, unless the algorithm needs to restart)
GATOC::evolve(fasta_file, vcf_file, :gen => gen, :pop_size => pop_size, :select_num => select_num, :c_mut => c_mut, :s_mut => s_mut,
 :save => save, :ran => ran, :loc => 'fragmented_genome_with_snps/arabidopsis_datasets',
 :start_pop => pop, :start_gen => restart_gen[0], :auc => 1, :auc_gen => 20, :restart_zero => restart_zero[0])