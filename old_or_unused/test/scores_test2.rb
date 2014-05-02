require '~/fragmented_genome_with_snps/lib/rearrangement_score.rb'
require '~/fragmented_genome_with_snps/lib/reform_ratio.rb'

fasta = ReformRatio::fasta_array("arabidopsis_datasets/ratio_dataset4/frags.fasta")

puts RearrangementScore::r_dist(fasta, fasta)
puts RearrangementScore::r_dist(fasta, fasta.shuffle)

puts RearrangementScore::gen_ham_dist(fasta, fasta)
puts RearrangementScore::gen_ham_dist(fasta, fasta.shuffle)
