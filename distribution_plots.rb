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

		ratios = WriteIt.file_to_floats_array("gen_#{i}_ratios.txt")
		SNPdist.plot_ratio(ratios, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", i, genome_length)

		hyp = SNPdist.hyp_snps(ratios, genome_length)
		SNPdist.plot_snps(hyp, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", i, genome_length, 'hyp', 
			'Approximated ratio of homozygous to heterozygous SNP density')

		hm = WriteIt.file_to_ints_array("gen_#{i}_hm.txt")
		SNPdist.plot_snps(hm, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", i, genome_length, 'hm',
			'Homozygous SNP density')

		ht = WriteIt.file_to_ints_array("gen_#{i}_ht.txt")
		SNPdist.plot_snps(ht, "fragmented_genome_with_snps/arabidopsis_datasets", "#{ARGV[0]}/#{ARGV[1]}", i, genome_length, 'ht',
			'Heterozygous SNP density')
	end
end