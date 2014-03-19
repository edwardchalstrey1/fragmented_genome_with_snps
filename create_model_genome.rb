#encoding: utf-8
require_relative 'lib/arabidopsis_c4_w_snps'
require_relative 'lib/write_it'
require_relative 'lib/snp_dist'

# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}"))

# Create the lists of homozygous and heterozygous SNPs
hm_str = SNPdist::hm
ht_str = SNPdist::ht
hm = hm_str[0]
ht = ht_str[0]
snp_pos = [ht, hm].flatten

puts "There are #{hm.length} homozygous SNPs"
puts "There are #{ht.length} heterozygous SNPs"

puts "Is there a SNP at the centre of the distribution? -- #{snp_pos.include?(10000000)}"

arabidopsis_c4 = ModelGenome::fasta_to_char_array("TAIR10_chr4.fasta")

contig_size = 10000 # 10-20kb
puts "Contig sizes = #{contig_size} to #{contig_size * 2}"
frags = ModelGenome::get_frags(arabidopsis_c4, contig_size)

puts "Arabidopsis chr4 length: #{arabidopsis_c4.length} bases"
puts "Fragmented seq   length: #{frags.join.length} = close enough? You decide."
puts "You have created #{frags.length} fragments of sizes #{contig_size}-#{contig_size*2}"

# Get the positions of the SNPs on fragments
pos_on_all = ModelGenome::pos_each_frag(snp_pos, frags)
pos_on_frags = pos_on_all[0]
snp_pos_all = pos_on_all[1]

fastaformat_array = ModelGenome::fasta_array(frags)
fastaformat_array_shuf = fastaformat_array.shuffle #shuffle it to show that the order doesn't need to be conserved when working out density later on

vcf = ModelGenome::vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)

WriteIt::write_data("arabidopsis_datasets/#{ARGV[0]}/frags.fasta", fastaformat_array)
WriteIt::write_data("arabidopsis_datasets/#{ARGV[0]}/snps.vcf", vcf)
WriteIt::write_data("arabidopsis_datasets/#{ARGV[0]}/frags_shuffled.fasta", fastaformat_array_shuf)
WriteIt::write_txt("arabidopsis_datasets/#{ARGV[0]}/info", [hm_str[1], ht_str[1], "Contig size = #{contig_size}"])