#encoding: utf-8
require "~/fragmented_genome_with_snps/lib/rearrangement_score"

## small_dataset2 ##
def get_perm(run, gen)
	perm = []
	IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/small_dataset2/run#{run}/Gen#{gen}/best_permutation.txt") { |line| perm << line.gsub(/\n/,'') }
	return perm[1..-1]
end
original = []; x = 1
53.times{original << "frag#{x}"; x+=1}

run2_gen3044 = get_perm(2, 3044)
run3_gen37 = get_perm(3, 37)
run4_gen10 = get_perm(4, 10)

puts RearrangementScore::dev_dist(original, run2_gen3044)
puts RearrangementScore::dev_dist(original, run3_gen37)
puts RearrangementScore::dev_dist(original, run4_gen10)
puts
puts RearrangementScore::dev_dist(run3_gen37, run2_gen3044)
puts RearrangementScore::dev_dist(run3_gen37, run4_gen10)
puts RearrangementScore::dev_dist(run2_gen3044, run4_gen10)
