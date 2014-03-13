#encoding: utf-8
require_relative 'lib/arabidopsis_c4_w_snps'
require_relative 'lib/write_it'

=begin
	
Use a uniform distribution for the heterozygous snps http://astrostatistics.psu.edu/su07/R/html/stats/html/Uniform.html
	
=end


# make the directory to put data files into
Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}"))

# Create the lists of homozygous and heterozygous SNPs
snpz = ModelGenome::snp_dist
snp_pos = snpz[0]
hm = snpz[1]
ht = snpz[2]

puts "Is there a SNP at the causative mutation position? -- #{snp_pos.include?(10000000)}"

arabidopsis_c4 = ModelGenome::fasta_to_char_array("TAIR10_chr4.fasta")

contig_size = 10000 # 10-20kb
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