#encoding: utf-8
require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'
require_relative 'lib/rearrangement_score'
require_relative 'lib/write_it'

if ARGV[0] == nil
	puts 'Wrong number of command line arguments!'
else
	vcf_file = "arabidopsis_datasets/#{ARGV[0]}/snps.vcf"
	ordered_fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags.fasta"
	ordered_fasta = ReformRatio::fasta_array(ordered_fasta_file)
	shuffled_population = GATOC::initial_population(ordered_fasta, 100) # number of scores to test

	scores = []
	shuffled_population.each do |perm| #permutation
		scores << RearrangementScore::rearrangement_score(ordered_fasta, perm)
	end
	worst_possible = RearrangementScore::rearrangement_score(ordered_fasta, ordered_fasta.reverse)
	average = scores.inject(:+)/scores.length
	best = scores.sort[0]
	worst = scores.sort[-1]
	perfect = RearrangementScore::rearrangement_score(ordered_fasta, ordered_fasta)
	WriteIt::write_txt("arabidopsis_datasets/#{ARGV[0]}/shuffled_ordinal_similarities", ["Perfect = #{perfect}", "Best = #{best}", "Average = #{average}", "Worst = #{worst}", "Worst Possible = #{worst_possible}", scores].flatten)
end
