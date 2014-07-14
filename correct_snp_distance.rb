require_relative 'lib/reform_ratio'
require_relative 'lib/GATOC'

fasta = ReformRatio.fasta_array('arabidopsis_datasets/10K_dataset4/frags.fasta')
snp_data = ReformRatio.get_snp_data('arabidopsis_datasets/10K_dataset4/snps.vcf')
ht, hm = ReformRatio.perm_pos(fasta, snp_data)

puts GATOC.snp_distance(hm)