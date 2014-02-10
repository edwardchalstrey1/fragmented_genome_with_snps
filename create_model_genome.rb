#encoding: utf-8
require 'rubygems'
require 'rinruby'
require 'bio-samtools'
require 'bio'
require_relative 'lib/arabidopsis_c4_w_snps'

Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/"+ARGV[0].to_s)) # make the directory to put data files into

snpz = ModelGenome::normal_dist
snp_pos = snpz[0]
hm = snpz[1]
ht = snpz[2]

puts "Is there a SNP at the causative mutation position? -- #{snp_pos.include?(10000000)}"

arabidopsis_c4 = ModelGenome::fasta_to_char_array("TAIR10_chr4.fasta")
snp_sequence = ModelGenome::snp_seq(arabidopsis_c4, snp_pos)

contig_size = 10000 # 10-20kb
frags = ModelGenome::get_frags(snp_sequence, contig_size)

puts "Arabidopsis chr4 length: #{arabidopsis_c4.length} bases"
puts "Fragmented seq   length: #{frags.join.length} = close enough? You decide."
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"

pos_on_all = ModelGenome::pos_each_frag(snp_pos, frags)
pos_on_frags = pos_on_all[0]
snp_pos_all = pos_on_all[1]

fastaformat_array = ModelGenome::fasta_array(frags)
fastaformat_array_shuf = fastaformat_array.shuffle #shuffle it to show that the order doesn't need to be conserved when working out density later on

vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

ModelGenome::write_data(fastaformat_array, 'frags.fasta', ARGV[0])
ModelGenome::write_data(vcf, 'snps.vcf', ARGV[0])
ModelGenome::write_data(fastaformat_array_shuf, 'frags_shuffled.fasta', ARGV[0])