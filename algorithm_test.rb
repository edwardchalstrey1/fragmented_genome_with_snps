#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'
require_relative 'lib/model_genome'
require_relative 'lib/snp_dist'
require_relative 'lib/write_it'

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

if dataset == nil || run == nil
	puts 'Wrong number of command line arguments! Specify a dataset and test-set'
else
	## Files ##
	vcf_file = "arabidopsis_datasets/#{dataset}/snps.vcf"
	fasta_file = "arabidopsis_datasets/#{dataset}/frags_shuffled.fasta"
	location = 'fragmented_genome_with_snps/arabidopsis_datasets'


	## Comparable ratio ##
	genome_length = ReformRatio::genome_length(fasta_file)
	hm = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/hm_snps.txt") # we can use the SNPs from the model genome to make example ratio
	ht = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/ht_snps.txt")
	fratio_breaks = SNPdist::fratio(hm, ht, div, genome_length) # frequency ratio array # 10,000 for 10K dataset
	comparable_ratio = SNPdist::hyp_snps(fratio_breaks, div, genome_length) # hypothetical snp positions array


	## Initial plots of correct ratio distribution ##
	Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{run}")) # make the directory to put data files into
	SNPdist::plot_hyp(comparable_ratio, location, "#{dataset}/#{run}", 'correct', genome_length)
	SNPdist::plot_dens(hm, ht, location, "#{dataset}/#{run}", genome_length, 'density_vector_ratio')


	## Plot of correctly ordered contigs ##
	correct_fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")
	snp_data = ReformRatio::get_snp_data(vcf_file)
	snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], correct_fasta) # array of no. of snps per frag in same order as correct_fasta
	pos_n_info = ReformRatio::get_positions(correct_fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) #get snp positions for each frag in array of arrays
	actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(correct_fasta)[1])
	het_snps, hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])
	fratio_breaks_ordered = SNPdist::fratio(hom_snps, het_snps, div, genome_length) # frequency ratio array # 10,000 for 10K dataset
	ratio_ordered = SNPdist::hyp_snps(fratio_breaks_ordered, div, genome_length) # hypothetical snp positions array
	SNPdist::plot_hyp(comparable_ratio, location, "#{dataset}/#{run}", 'correct_ordered_contigs', genome_length)
	SNPdist::plot_dens(hm, ht, location, "#{dataset}/#{run}", genome_length, 'density_vector_ratio_ordered_contigs')


	## Average fitness of correctly ordered contigs ##
	average_fitness_correct = []
	10.times do
		average_fitness_correct << GATOC::fitness(correct_fasta, snp_data, 'x', comparable_ratio, 'loc', 'dataset', 'run', div, genome_length)[0]
	end
	average = (average_fitness_correct.inject(:+))/10.0
	WriteIt::write_txt("arabidopsis_datasets/#{dataset}/#{run}/average_fitness_correct", [average])


	## Run the algorithm ##
	GATOC::evolve(fasta_file, vcf_file, :gen => gen, :pop_size => pop_size, :select_num => select_num, :c_mut => c_mut, :s_mut => s_mut,
	 :save => save, :ran => ran, :loc => 'fragmented_genome_with_snps/arabidopsis_datasets', :comparable_ratio => comparable_ratio, 
	 :div => div, :genome_length => genome_length, :end => average, :gen_end => 10)
end