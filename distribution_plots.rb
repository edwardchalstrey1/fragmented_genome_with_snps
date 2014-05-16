#encoding: utf-8
require_relative 'lib/score_plots/example_perms'
require_relative 'lib/snp_dist'
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'

=begin
	
ARGV's: 0 = dataset
		1 = run
		2 = generation (number)
=end

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")

Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{ARGV[2]}_lists")) do
	perm_ratio = WriteIt.file_to_floats_array("gen_#{ARGV[2]}_hyp.txt")
	SNPdist::plot_hyp(perm_ratio, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", ARGV[2], genome_length)
end
