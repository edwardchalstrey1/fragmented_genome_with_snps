#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'
require_relative 'lib/fitness_score'
require_relative 'lib/snp_dist'
require_relative 'lib/write_it'
require_relative 'lib/score_plots/score_plots'
require_relative 'lib/score_plots/example_perms'

dataset = ARGV[0]
run = ARGV[1]
gen =  ARGV[2].to_i # Number of generations
pop_size = ARGV[3].to_i # Population size
select_num = ARGV[4].to_i # Number of permutations to select from each generation
c_mut = ARGV[5].to_i # Number of chunk mutants (see GATOC and PMeth gem) in each new population
s_mut = ARGV[6].to_i # Number of swap mutants (see GATOC and PMeth gem) in each new population
save = ARGV[7].to_i # Number of permutations to save from each generation
ran = ARGV[8].to_i # Number of random permutations in each generation
div = ARGV[9].to_f # Number of divisions in the genome, at which to calculate the SNP frequencies
restart = ARGV[10] # Tells the algorithm to continue from the most recent generation if it has stopped

## Files ##
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"
location = 'fragmented_genome_with_snps/arabidopsis_datasets'

snp_data = ReformRatio::get_snp_data(vcf_file)
fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta") # correct permutation

### Directories ###

Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) # make the directory to put permutation files into
Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}/Gencorrect_lists")) # make the directory to put correct permutation files into

## Comparable ratio ## TODO a comparable ratio that doesn't use the known distributions
genome_length = ReformRatio::genome_length(fasta_file)
snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) # array of no. of snps per frag in same order as fasta
pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) # get snp positions for each frag in array of arrays
actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
ht, hm = ReformRatio::het_hom(actual_pos, pos_n_info[1])
WriteIt::write_txt("arabidopsis_datasets/#{dataset}/hm_snps", hm)
WriteIt::write_txt("arabidopsis_datasets/#{dataset}/ht_snps", ht)
hom_count = FitnessScore::count(hm, div, genome_length)
het_count = FitnessScore::count(ht, div, genome_length)
comparable_ratio = FitnessScore::ratio(hom_count, het_count)

SNPdist.plot_ratio(comparable_ratio, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length)

hyp = SNPdist.hyp_snps(comparable_ratio, genome_length)
SNPdist.plot_snps(hyp, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, "hyp_#{div/1000}Kdiv", 
	'Approximated ratio of homozygous to heterozygous SNP density')

SNPdist.plot_snps(hm, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, "hm_#{div/1000}Kdiv",
	'Homozygous SNP density')

SNPdist.plot_snps(ht, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, "ht_#{div/1000}Kdiv",
	'Heterozygous SNP density')

pop, restart_gen = [], []
if restart == nil
	fitness_correct = GATOC::fitness(fasta, snp_data, comparable_ratio, div, genome_length)[0]
	WriteIt::write_txt("arabidopsis_datasets/#{dataset}/#{run}/fitness_correct", [fitness_correct])
	restart_gen << 0
	pop = nil
else ## For restarts ##
	Dir.chdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) do # Selecting the generation to restart from
		dirs = Dir.glob('*').select {|f| File.directory? f}
		dirs.delete("Gencorrect_lists")
		restart_gen << (dirs.length/2)-1
	end
	id_pop = MetricPlot::get_perms(1, restart_gen[0], 0, dataset, run).flatten(1) # should be population array of permutations (arrays of ids)
	id_pop.each do |ids|
		pop << [ExamplePerms::fasta_p_id(fasta, ids), 'type']
	end
end

# Run the algorithm ##
GATOC::evolve(fasta_file, vcf_file, :gen => gen, :pop_size => pop_size, :select_num => select_num, :c_mut => c_mut, :s_mut => s_mut,
 :save => save, :ran => ran, :loc => 'fragmented_genome_with_snps/arabidopsis_datasets', :comparable_ratio => comparable_ratio, 
 :div => div, :genome_length => genome_length, :start_pop => pop, :start_gen => restart_gen[0], :auc => 1, :auc_gen => 20)