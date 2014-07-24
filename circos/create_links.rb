require_relative '../lib/write_it'
require_relative '../lib/reform_ratio'

class CircosLinks

	def self.original_order
		original_fasta = ReformRatio::fasta_array("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/frags.fasta")
		return ReformRatio::fasta_id_n_lengths(original_fasta)[0]
	end

	def self.best_perms(start, last, step)
		range = []
		(start..last).step(step) do |n|
			range << n
		end
		best_perms = [] # best permutation for each generation in ascending order
		range.each do |i|
			permutation = []
			IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{i}/best_permutation.txt") { |line| permutation << line.gsub(/\n/,'') }
			best_perms << permutation[1..-1]
		end
		return best_perms, range
	end

	def self.links(best_perms, original_order)
		links = []
		chrom = 2
		best_perms.each do |perm|
			f = 1 # position of the fragment in this permutation
			perm.each do |frag|
				pos = original_order.index(frag) + 1 # position of frag in original order
				links << "p1 #{pos} #{pos+1} p#{chrom} #{f} #{f+1} score=#{pos}" # colour score is based on position in original order, p for permutation
				f+=1
			end
			chrom+=1
		end
		links
	end

	def self.karyotype(best_perms, range)
		frag_num = best_perms[0].length
		karyotype = ["chr - p1 Correct 0 #{frag_num} blue"]
		x = 0
		best_perms.each do |perm|
			karyotype << "chr - p#{x+2} Gen_#{range[x]} 0 #{frag_num} green"
			x+=1
		end
		karyotype
	end

	def self.config(start, last, frag_num)
		config = ["karyotype = data/karyotype_#{ARGV[0]}_#{ARGV[1]}_#{start}-#{last}.txt",
		'chromosomes_radius  = hs1:0.95r',
		'<links>',
		'<link>',
		"file = data/links_#{ARGV[0]}_#{ARGV[1]}_#{start}-#{last}.txt",
		'radius        = 0.99r',
		'bezier_radius = 0.1r',
		'color         = blue',
		'thickness     = 1',
		'<rules>',
		'<rule>', # for all links
		'condition = 1',
		"color = eval( sprintf('hue%03d',remap_int(var(score),1,#{frag_num},0,060)))", # frag_num is the number of fragments
		'</rule>',
		'</rules>',
		'</link>',
		'</links>',
		'<<include ideogram.conf>>',
		'<image>',
		'<<include etc/image.conf>>',                
		'</image>',
		'<<include etc/colors_fonts_patterns.conf>>',
		'<<include etc/housekeeping.conf>> ']
		return config
	end

	def self.make_links(start, last, step)
		original = original_order
		best_n_range = best_perms(start, last, step)
		best = best_n_range[0]
		range = best_n_range[1]
		WriteIt::write_txt("circos/data/links_#{ARGV[0]}_#{ARGV[1]}_#{start}-#{last}", links(best, original))
		WriteIt::write_txt("circos/data/karyotype_#{ARGV[0]}_#{ARGV[1]}_#{start}-#{last}", karyotype(best, range))
		WriteIt::write_data("circos/#{ARGV[0]}_#{ARGV[1]}_#{start}-#{last}.conf", config(start, last, original.length))
	end
end

# ARGV[0] == dataset
# ARGV[1] == run
# ARGV[2] == lower end of range of generations for circos plot
# ARGV[3] == upper end
# ARGV[4] == step of range