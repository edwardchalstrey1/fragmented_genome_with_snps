#encoding: utf-8
require_relative 'lib/score_plots/example_perms.rb'
require_relative 'lib/score_plots/score_plots.rb'
require_relative 'lib/snp_dist'

### Example plots ###

location = "fragmented_genome_with_snps/arabidopsis_datasets"
dataset = ARGV[0]

div = ARGV[1].to_f # Number of divisions at which to calculate SNP density in permutation
n = ARGV[2].to_i # Number of permutations in each population (mutants of the original order, and random shuffles)

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{dataset}/snps.vcf")
fasta = ReformRatio::fasta_array("arabidopsis_datasets/#{dataset}/frags.fasta")

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{dataset}/frags.fasta")
hm = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/hm_snps.txt") # we can use the SNPs from the model genome to make example ratio
ht = WriteIt::file_to_ints_array("#{Dir.home}/#{location}/#{dataset}/ht_snps.txt")

hom_count = FitnessScore::count(hm, div, genome_length)
het_count = FitnessScore::count(ht, div, genome_length)
comparable_ratio = FitnessScore::ratio(hom_count, het_count)

example_perms = ExamplePerms::get_perms(fasta, n, snp_data, comparable_ratio, div, genome_length) 

Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{div/1000}K_mutant_perm_examples"))

MetricPlot::metric_plot(0, 2, "dev", "#{n}_chunk_swap_shuf_dev", example_perms, dataset, "#{div/1000}K_mutant_perm_examples")
MetricPlot::metric_plot(0, 2, "sq", "#{n}_chunk_swap_shuf_sq", example_perms, dataset, "#{div/1000}K_mutant_perm_examples")
MetricPlot::metric_plot(0, 2, "ham", "#{n}_chunk_swap_shuf_ham", example_perms, dataset, "#{div/1000}K_mutant_perm_examples")
MetricPlot::metric_plot(0, 2, "r", "#{n}_chunk_swap_shuf_r", example_perms, dataset, "#{div/1000}K_mutant_perm_examples")
MetricPlot::metric_plot(0, 2, "lcs", "#{n}_chunk_swap_shuf_lcs", example_perms, dataset, "#{div/1000}K_mutant_perm_examples")
MetricPlot::metric_plot(0, 2, "kt", "#{n}_chunk_swap_shuf_kt", example_perms, dataset, "#{div/1000}K_mutant_perm_examples")

dirs = ['chunk_mutants', 'swap_mutants', 'shuffled']
x = 0
example_perms.each do |population|
	y = 1
	Dir.mkdir(File.join(Dir.home, "#{location}/#{dataset}/#{div/1000}K_mutant_perm_examples/#{dirs[x]}"))
	population.each do |permutation|
		WriteIt::write_txt("#{Dir.home}/#{location}/#{dataset}/#{div/1000}K_mutant_perm_examples/#{dirs[x]}/permutation#{y}", permutation)
		fasta_perm = ExamplePerms::fasta_p_id(fasta, permutation)
		correlation, hom_snps, het_snps, perm_ratio = GATOC::fitness(fasta_perm, snp_data, comparable_ratio, div, genome_length)
		SNPdist::plot_ratio2(perm_ratio, location, "#{dataset}/#{div/1000}K_mutant_perm_examples", dirs[x], genome_length, "permutation#{y}")

		hyp = SNPdist.hyp_snps(perm_ratio, genome_length)
		SNPdist.plot_snps2(hyp, "fragmented_genome_with_snps/arabidopsis_datasets", "#{dataset}/#{div/1000}K_mutant_perm_examples",
			dirs[x], genome_length, "permutation#{y}_hyp", 
			'Approximated ratio of homozygous to heterozygous SNP density')
		y+=1
	end
	x+=1
end
