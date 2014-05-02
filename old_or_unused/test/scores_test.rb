require_relative 'lib/rearrangement_score.rb'
require 'pp'
require_relative 'lib/write_it'

orig = (1..4).to_a
a = []
orig.permutation.map(&:join).each {|i| a << i.split(//)}
x = 0
a.each do |i|
	a[x] = []
	i.each {|j| a[x] << j.to_i}
	x+=1
end

=begin
	
a is all the possible permutations of 1,2,3,4
	
=end

dev, sq_dev, ham, mod_ham, r, lcs, kt = [], [], [], [], [], [], []
a.each do |perm|
	dev << RearrangementScore::dev_dist(orig, perm)
	sq_dev << RearrangementScore::sq_dev_dist(orig, perm)
	ham << RearrangementScore::gen_ham_dist(orig, perm)
	mod_ham << RearrangementScore::mod_ham_dist(orig, perm)
	r << RearrangementScore::r_dist(orig, perm)
	lcs << RearrangementScore::lcs(orig, perm)
	kt << RearrangementScore::kendalls_tau(orig, perm)
end

scores = dev.zip(sq_dev).zip(ham).zip(mod_ham).zip(r).zip(lcs).zip(kt).flatten
perms = []
a.each {|perm| perms << perm.join.to_s}

WriteIt::write_txt('test/1-4_metric_scores', scores)
WriteIt::write_txt('test/1-4_perms', perms)
#WriteIt::write_txt('test/1-4_perms_arrays', a)
