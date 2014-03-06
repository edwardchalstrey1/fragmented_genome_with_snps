#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/write_it'
require '~/fragmented_genome_with_snps/lib/rearrangement_score'
require 'rinruby'

def original_order
	original_order = []
	IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/correct_permutation.txt") { |line| original_order << line }
	original_order[1..-1]
end

def get_perms(pop_size)
	all_perms = []
	pop_size.times do
		x = 0
		pop = []
		Dir.entries("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{x}").each do |ptxt|
			if ptxt.include? '.txt'
				perm = []
				IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{x}/#{ptxt}") { |line| perm << line }
				pop << perm
			end
		end
		all_perms << pop
		x+=1
	end
	return all_perms
end

def get_metrics(all_perms)
	orig = original_order
	gen_pop_met = []
	all_perms.each do |pop|
		pop_met = []
		pop.each do |perm|
			metrics = {}
			metrics[:fit] = perm[0]
			metrics[:dev] = RearrangementScore::dev_dist(orig, perm[1..-1])
			metrics[:sq] = RearrangementScore::sq_dev_dist(orig, perm[1..-1])
			metrics[:ham] = RearrangementScore::gen_ham_dist(orig, perm[1..-1])
			metrics[:mod] = RearrangementScore::mod_ham_dist(orig, perm[1..-1])
			metrics[:r] = RearrangementScore::r_dist(orig, perm[1..-1])
			metrics[:lcs] = RearrangementScore::lcs(orig, perm[1..-1])
			metrics[:kt] = RearrangementScore::kendalls_tau(orig, perm[1..-1])
			pop_met << metrics
		end
		gen_pop_met << pop_met
	end
	return gen_pop_met
end

all_perms = get_perms(29) # [][][0] is the fitness of the permutation [][][1..-1] is the permutation (frag ids)
gen_pop_met = get_metrics(all_perms)

puts gen_pop_met[0][0]

