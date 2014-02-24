require_relative 'lib/write_it'

def original_order
	original_order = []
	IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/correct_permutation.txt") { |line| original_order << line }
	original_order
end


def best_perms
	best_perms = [] # best permutation for each generation in ascending order
	y = 100
	x = 1
	y.times do
		permutation = []
		IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{x}_best_permutation.txt") { |line| permutation << line } # GEN NEEDS TO CHANGE
		best_perms << permutation
		x+=1
	end
	best_perms
end

def links (best_perms, original_order, lrange, urange)
	links = []
	chrom = 2
	best_perms[(lrange-1)..(urange-1)].each do |perm|
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

def karyotype (best_perms, original_order, lrange, urange)
	frag_num = best_perms[0].length
	karyotype = ["chr - p1 Correct 0 #{frag_num} chr1"]
	x = 2
	best_perms[(lrange-1)..(urange-1)].each do |perm|
		karyotype << "chr - p#{x} Gen_#{x-1} 0 #{frag_num} green"
		x+=1
	end
	karyotype
end

original = original_order
best = best_perms

WriteIt::write_txt("circos/data/links_#{ARGV[0]}_#{ARGV[1]}_#{ARGV[2]}-#{ARGV[3]}", links(best, original, ARGV[2].to_i, ARGV[3].to_i))
WriteIt::write_txt("circos/data/karyotype_#{ARGV[0]}_#{ARGV[1]}_#{ARGV[2]}-#{ARGV[3]}", karyotype(best, original, ARGV[2].to_i, ARGV[3].to_i))