#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/snp_dist'
require_relative 'lib/fitness_score'

=begin
	
ARGV's: 0 = dataset
		1 = run
		2 = generation (number)
		3 = div
=end
dataset = ARGV[0]
run = ARGV[1]
gen = ARGV[2]
div = ARGV[3].to_f

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")

Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}/#{run}/Gencorrect_lists")) # make the directory to put correct permutation files into

Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}")) do
	hm = WriteIt.file_to_ints_array("hm_snps.txt")
	ylim_hm = SNPdist.get_ylim(hm, genome_length, 'density')
	ht = WriteIt.file_to_ints_array("ht_snps.txt")
	ylim_ht = SNPdist.get_ylim(ht, genome_length, 'density')
	hom_count = FitnessScore::count(hm, div, genome_length)
	het_count = FitnessScore::count(ht, div, genome_length)
	puts hom_count # TODO what is going on here?
	puts het_count
	comparable_ratio = FitnessScore::ratio(hom_count, het_count)
	ylim_ratio = SNPdist.get_ylim(comparable_ratio, genome_length, 'ratio')
	SNPdist.plot_ratio(comparable_ratio, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, ylim_ratio)

	hyp = SNPdist.hyp_snps(comparable_ratio, genome_length)
	ylim_hyp = SNPdist.get_ylim(hyp, genome_length, 'density')
	SNPdist.plot_snps(hyp, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, "hyp_#{div/1000}Kdiv", 
		'Approximated ratio of homozygous to heterozygous SNP density', ylim_hyp)

	SNPdist.plot_snps(hm, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, "hm_#{div/1000}Kdiv",
		'Homozygous SNP density', ylim_hm)

	SNPdist.plot_snps(ht, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", 'correct', genome_length, "ht_#{div/1000}Kdiv",
		'Heterozygous SNP density', ylim_ht)
end

Array(0..gen.to_i).each do |i|
	Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{dataset}/#{run}/Gen#{i}_lists")) do

		ratios = WriteIt.file_to_floats_array("gen_#{i}_ratios.txt")
		SNPdist.plot_ratio(ratios, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", i, genome_length, ylim_ratio)

		hyp = SNPdist.hyp_snps(ratios, genome_length)
		SNPdist.plot_snps(hyp, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", i, genome_length, 'hyp', 
			'Approximated ratio of homozygous to heterozygous SNP density', ylim_hyp)

		hm = WriteIt.file_to_ints_array("gen_#{i}_hm.txt")
		SNPdist.plot_snps(hm, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", i, genome_length, 'hm',
			'Homozygous SNP density', ylim_hm)

		ht = WriteIt.file_to_ints_array("gen_#{i}_ht.txt")
		SNPdist.plot_snps(ht, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{run}", i, genome_length, 'ht',
			'Heterozygous SNP density', ylim_ht)
	end
end

