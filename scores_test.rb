require_relative 'lib/rearrangement_score.rb'
require 'pp'
require_relative 'lib/write_it'

vp = [1,2,3]
vo = [1,2,3]

#puts RearrangementScore::rearrangement_score(vo, vp)
#puts RearrangementScore::metric_2(vo, vp)

orig = (1..4).to_a
a = []
orig.permutation.map(&:join).each {|i| a << i.split(//)}
x = 0
a.each do |i|
	a[x] = []
	i.each {|j| a[x] << j.to_i}
	x+=1
end

m1_scores, m2a_scores, m2b_scores, m2c_scores, m12av_scores = [], [], [], [], []
a.each do |perm|
	m1_scores << RearrangementScore::rearrangement_score(orig, perm)
	m2a_scores << RearrangementScore::metric_2(orig, perm, 'a')
	m2b_scores << RearrangementScore::metric_2(orig, perm, 'b')
	m2c_scores << RearrangementScore::metric_2(orig, perm, 'c')
	m12av_scores << RearrangementScore::metric_1_2_av(orig, perm)
end

scores = m1_scores.zip(m2a_scores).zip(m2b_scores).zip(m2c_scores).zip(m12av_scores).flatten

perms = []
a.each {|perm| perms << perm.join.to_s}

WriteIt::write_txt('test/1-4_metric_scores', scores)
WriteIt::write_txt('test/1-4_perms', perms)
#WriteIt::write_txt('test/1-4_perms_arrays', a)
