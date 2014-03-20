#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/snp_dist'
require '~/fragmented_genome_with_snps/lib/reform_ratio'
require '~/fragmented_genome_with_snps/lib/GATOC'

vcf_file = "arabidopsis_datasets/10K_dataset1/snps.vcf"
snp_data = ReformRatio::get_snp_data(vcf_file)
fasta_file = "arabidopsis_datasets/10K_dataset1/frags.fasta"
fasta = ReformRatio::fasta_array(fasta_file)

fasta_shuf = fasta.shuffle

def hyp_snps(fasta, snp_data)
	snps_per_frag = ReformRatio::snps_per_fasta_frag(snp_data[2], fasta) #array of no. of snps per frag in same order as fasta
	pos_n_info = ReformRatio::get_positions(fasta, snp_data[0], snp_data[1], snps_per_frag, snp_data[3]) #get snp positions for each frag in array of arrays
	actual_pos = ReformRatio::total_pos(pos_n_info[0], ReformRatio::fasta_id_n_lengths(fasta)[1])
	het_hom_snps = ReformRatio::het_hom(actual_pos, pos_n_info[1])
	fratio_breaks_real = SNPdist::fratio(het_hom_snps[1], het_hom_snps[0], 10000)
	return SNPdist::hyp_snps(fratio_breaks_real, 10000), het_hom_snps[0], het_hom_snps[1]
end

# hm = SNPdist::hm[0]
# ht = SNPdist::ht[0]

# fratio_breaks = SNPdist::fratio(hm, ht, 10000)
# hyp = SNPdist::hyp_snps(fratio_breaks, 10000)

# hyp_orig = hyp_snps(fasta, snp_data)
# hyp_shuf = hyp_snps(fasta_shuf, snp_data)

# SNPdist::plot_hyp(hyp_shuf[0])
# SNPdist::plot_dens(hyp_shuf[2], hyp_shuf[1])
# puts SNPdist::qq_cor(hyp, hyp_shuf[0])

# SNPdist::plot_hyp(hyp_orig[0])
# SNPdist::plot_dens(hyp_orig[2], hyp_orig[1])
# puts SNPdist::qq_cor(hyp, hyp_orig[0])



### Look at these Scores ###

puts GATOC::fitness(fasta, snp_data, 'same')
puts GATOC::fitness(fasta, snp_data, 'diff')

puts GATOC::fitness(fasta_shuf, snp_data, 'same')
puts GATOC::fitness(fasta_shuf, snp_data, 'diff')