require_relative 'lib/rearrangement_score.rb'
require 'pp'

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

m1_scores, m2_scores = [], []
a.each do |perm|
	m1_scores << RearrangementScore::rearrangement_score(orig, perm)
	m2_scores << RearrangementScore::metric_2(orig, perm)
end

puts m1_scores
puts
puts m2_scores # now write these to txts or something
