#encoding: utf-8
require_relative 'lib/reform_ratio'

vcf_file = "arabidopsis_datasets/#{ARGV[0]}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags_shuffled.fasta"

ordered_fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags.fasta"
ordered_fasta = ReformRatio::fasta_array(ordered_fasta_file)

# fasta_file, vcf_file, gen, pop_size, select_num, mut_num, save, ran, ordered_fasta, figures
GATOC::evolve(fasta_file, vcf_file, 10, 20, 10, 2, 1, 1, ordered_fasta, "no figures")
