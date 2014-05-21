#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/write_it'
require_relative 'lib/snp_dist'

=begin
	
ARGV's: 0 = dataset
		1 = run
		2 = generation (number)
=end

genome_length = ReformRatio::genome_length("arabidopsis_datasets/#{ARGV[0]}/frags.fasta")

Array(0..ARGV[2].to_i).each do |i|
	Dir.chdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{i}_lists")) do
		perm_ratio = WriteIt.file_to_floats_array("gen_#{i}_hyp.txt")
		SNPdist::plot_hyp(perm_ratio, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", i, genome_length)
	end
end