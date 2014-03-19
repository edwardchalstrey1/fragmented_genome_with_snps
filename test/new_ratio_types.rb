#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/snp_dist'
require '~/fragmented_genome_with_snps/lib/reform_ratio'
require '~/fragmented_genome_with_snps/lib/GATOC'

vcf_file = "arabidopsis_datasets/10K_dataset1/snps.vcf"
snp_data = ReformRatio::get_snp_data(vcf_file)
fasta_file = "arabidopsis_datasets/10K_dataset1/frags_shuffled.fasta"
fasta = ReformRatio::fasta_array(fasta_file)

# hm = SNPdist::hm[0]
# ht = SNPdist::ht[0]

# fratio_breaks = SNPdist::fratio(hm, ht, 10000)
# hyp = SNPdist::hyp_snps(fratio_breaks, 10000)

# SNPdist::plot_hyp(hyp,hm,ht)

# sample_ratio = SNPdist::sample_ratio(hm,ht)

### Look at these Scores ###

puts GATOC::fitness(fasta, snp_data, 'same')
puts GATOC::fitness(fasta, snp_data, 'diff')

puts GATOC::fitness(fasta.shuffle, snp_data, 'same')
puts GATOC::fitness(fasta.shuffle, snp_data, 'diff')