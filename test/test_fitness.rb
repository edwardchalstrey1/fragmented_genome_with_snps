require '~/fragmented_genome_with_snps/lib/reform_ratio.rb'
require '~/fragmented_genome_with_snps/lib/GATOC.rb'


snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf")

x = 1
["arabidopsis_datasets/#{ARGV[0]}/frags_shuffled.fasta", "arabidopsis_datasets/#{ARGV[0]}/frags.fasta"].each do |fasta_file|
	fasta = ReformRatio::fasta_array(fasta_file)
	puts "#{x} cor: #{GATOC::fitness(fasta, snp_data, 'same')}"
	puts "#{x} kol: #{GATOC::fitness(fasta, snp_data, 'kol')}"
	x+=1
end