require_relative 'lib/write_it'

def original_order
	original_order = []
	IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/correct_permutation.txt") { |line| original_order << line }
	original_order
end


def best_perms(lrange, urange, step)
	range = []
	(lrange..urange).step(step) do |n|
		range << n
	end
	best_perms = [] # best permutation for each generation in ascending order
	range.each do |i|
		permutation = []
		IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{i}_best_permutation.txt") { |line| permutation << line }
		# permutation = permutation[1..-1]
		best_perms << permutation 
	end
	return best_perms, range
end

def links(best_perms, original_order)
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

def karyotype(best_perms, range)
	frag_num = best_perms[0].length
	karyotype = ["chr - p1 Correct 0 #{frag_num} blue"]
	x = 0
	best_perms.each do |perm|
		karyotype << "chr - p#{x+2} Gen_#{range[x]} 0 #{frag_num} green"
		x+=1
	end
	karyotype
end

original = original_order
best_n_range = best_perms(ARGV[2].to_i, ARGV[3].to_i, ARGV[4].to_i)
best = best_n_range[0]
range = best_n_range[1]

WriteIt::write_txt("circos/data/links_#{ARGV[0]}_#{ARGV[1]}_#{ARGV[2]}-#{ARGV[3]}", links(best, original))
WriteIt::write_txt("circos/data/karyotype_#{ARGV[0]}_#{ARGV[1]}_#{ARGV[2]}-#{ARGV[3]}", karyotype(best, range))

# ARGV[0] == dataset
# ARGV[1] == run
# ARGV[2] == lower end of range of generations for circos plot
# ARGV[3] == upper end
# ARGV[4] == step of range