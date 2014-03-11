require '~/fragmented_genome_with_snps/lib/reform_ratio.rb'
require '~/fragmented_genome_with_snps/lib/GATOC.rb'

snp_data = ReformRatio::get_snp_data("arabidopsis_datasets/ratio_dataset4/snps.vcf")

fasta_s = ReformRatio::fasta_array("arabidopsis_datasets/ratio_dataset4/frags_shuffled.fasta")
fasta = ReformRatio::fasta_array("arabidopsis_datasets/ratio_dataset4/frags.fasta")
mut = GATOC::mutate(fasta)

def fasta_p(fasta, perm)
	perm_ids = []
	IO.foreach("arabidopsis_datasets/ratio_dataset4/run6/Gen#{perm}/best_permutation.txt") { |line| perm_ids << line.gsub(/\n/,'') }
	perm_ids = perm_ids[2..-1]
	ids = ReformRatio::fasta_id_n_lengths(fasta)[0]
	fasta_perm = []
	perm_ids.each do |id|
		fasta.each do |frag|
			if id == frag.entry_id
				fasta_perm << frag
			end
		end
	end
	fasta_perm
end

f212 = fasta_p(fasta, 212)

fastas = [fasta, fasta_s, mut, f212]

x = 1
fastas.each do |fasta|
	puts "#{x} cor: #{GATOC::fitness(fasta, snp_data, 'same')}"
	puts "#{x} kol: #{GATOC::fitness(fasta, snp_data, 'kol')}"
	puts "#{x} sp: #{GATOC::fitness(fasta, snp_data, 'spearman')}"
	x+=1
end