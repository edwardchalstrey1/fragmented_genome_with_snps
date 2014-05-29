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
gen =  ARGV[2].to_i
pop_size = ARGV[3].to_i
select_num = ARGV[4].to_i
c_mut = ARGV[5].to_i
s_mut = ARGV[6].to_i
save = ARGV[7].to_i
ran = ARGV[8].to_i
div = ARGV[9].to_f
restart = ARGV[10]

## Files ##
vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"
location = 'fragmented_genome_with_snps/arabidopsis_datasets'

snp_data = ReformRatio::get_snp_data(vcf_file)
correct_fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")

### Directories ###

Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) # make the directory to put permutation files into
Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}/Gencorrect_lists")) # make the directory to put correct permutation files into

## Comparable ratio ## TODO a comparable ratio that doesn't use the known distributions
genome_length = ReformRatio::genome_length(fasta_file)
hm = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/hm_snps.txt") # we can use the SNPs from the model genome to make example ratio
ht = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/ht_snps.txt")
hom_count = FitnessScore::count(hm, div, genome_length)
het_count = FitnessScore::count(ht, div, genome_length)
comparable_ratio = FitnessScore::ratio(hom_count, het_count) # TODO plot of correctly ordered contigs

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
	## Average fitness of correctly ordered contigs ##
	average_fitness_correct = []
	10.times do
		average_fitness_correct << GATOC::fitness(correct_fasta, snp_data, comparable_ratio, div, genome_length)[0]
	end
	average = (average_fitness_correct.inject(:+))/10.0
	WriteIt::write_txt("arabidopsis_datasets/#{dataset}/#{run}/average_fitness_correct", [average])

	restart_gen << 0
	pop = nil

else ## For restarts ##
	Dir.chdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) do # Selecting the generation to restart from
		dirs = Dir.glob('*').select {|f| File.directory? f}
		restart_gen << (dirs.length/2)-1
	end

	id_pop = MetricPlot::get_perms(1, restart_gen[0], 0, dataset, run).flatten(1) # should be population array of permutations (arrays of ids)
	id_pop.each do |ids|
		pop << [ExamplePerms::fasta_p_id(correct_fasta, ids), 'type']
	end
end

# Run the algorithm ##
GATOC::evolve(fasta_file, vcf_file, :gen => gen, :pop_size => pop_size, :select_num => select_num, :c_mut => c_mut, :s_mut => s_mut,
 :save => save, :ran => ran, :loc => 'fragmented_genome_with_snps/arabidopsis_datasets', :comparable_ratio => comparable_ratio, 
 :div => div, :genome_length => genome_length, :start_pop => pop, :start_gen => restart_gen[0])