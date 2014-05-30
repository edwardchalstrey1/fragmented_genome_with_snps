#encoding: utf-8
require_relative 'lib/GATOC'
require_relative 'lib/reform_ratio'
require_relative 'lib/snp_dist'
require_relative 'lib/quit_if'
require_relative 'lib/write_it'

ordered_fasta_array = ReformRatio.fasta_array('arabidopsis_datasets/10K_dataset3/frags.fasta')
genome_length = ReformRatio.genome_length('arabidopsis_datasets/10K_dataset3/frags.fasta')
snp_data = ReformRatio.get_snp_data('arabidopsis_datasets/10K_dataset3/snps.vcf')

puts 'Start'

def fits(fasta_array, snp_data, div, genome_length)

	snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta_array) #array of no. of snps per frag in same order as fasta
	pos_n_info = ReformRatio::get_positions(fasta_array, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) #get snp positions for each frag in array of arrays
	actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta_array)[1])
	het_snps, hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])
	fratio_breaks_perm = SNPdist::fratio(hom_snps, het_snps, div, genome_length) # frequency ratio array
	comparable_ratio = SNPdist::hyp_snps(fratio_breaks_perm, div, genome_length) # hypothetical snp positions array
	puts 'Got comparable_ratio'

	fitness_scores = []; x = 1
	10.times do
		fitness_scores << GATOC.fitness(fasta_array, snp_data, comparable_ratio, div, genome_length)[0]
		puts "got fitness #{x}"; x+=1
	end

	return fitness_scores.sort
end

div_1K_fits = fits(ordered_fasta_array, snp_data, 1000.0, genome_length)
div_10K_fits = fits(ordered_fasta_array, snp_data, 10000.0, genome_length)
div_100K_fits = fits(ordered_fasta_array, snp_data, 100000.0, genome_length)

slopes = []; x = 1
[div_1K_fits, div_10K_fits, div_100K_fits].each do |fits|
	slopes << QuitIf.quit(fits)
	puts "Got slope #{x}"; x+=1
end

WriteIt.write_txt('arabidopsis_datasets/10K_dataset3/slopes_div-1K-10K-100K', slopes)
puts 'It is written...'