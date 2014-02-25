require_relative 'lib/write_it'
require_relative 'lib/rearrangement_score'
require 'rinruby'
### Ordinal Similarity Through Generations ###

def original_order
	original_order = []
	IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/correct_permutation.txt") { |line| original_order << line }
	original_order
end

def best_perms(lrange, urange)
	range = *(lrange..urange)
	best_perms = [] # best permutation for each generation in ascending order
	range.each do |i|
		permutation = []
		IO.foreach("arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{i}_best_permutation.txt") { |line| permutation << line }
		# permutation = permutation[1..-1]
		best_perms << permutation 
	end
	return best_perms
end

def ord_sim_fig(best_perms, original_order, dataset, run)
	scores = []
	best_perms.each do |perm|
		scores << RearrangementScore::rearrangement_score(original_order, perm)
	end
	generations = *(1..best_perms.length)
	r_sesh = RinRuby.new(echo = false)
	r_sesh.assign 'generations', generations
	r_sesh.assign 'dataset', dataset
	r_sesh.assign 'run', run
	r_sesh.assign 'scores', scores
	r_sesh.eval 'png(paste("~/fragmented_genome_with_snps/arabidopsis_datasets/", dataset, "/", run, "/ord_sim_over_generations.png", sep=""))'
	r_sesh.eval 'plot(generations, scores, main="Line represents average score for shuffled permutations", xlab="Generation", ylab="Ordinal Similarity Score")'
	r_sesh.eval 'abline(512351, 0)' # needs to be average of randoms #TODO
	r_sesh.eval 'dev.off()'
	r_sesh.quit
end

original = original_order
best = best_perms(ARGV[2].to_i, ARGV[3].to_i)

ord_sim_fig(best, original, ARGV[0], ARGV[1])