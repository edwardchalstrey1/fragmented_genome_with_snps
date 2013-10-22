require 'munkres'
cost_matrix = [[3,2,1,0],[2,3,2,1],[1,2,3,2],[0,1,2,3]]
m = Munkres.new(cost_matrix)
p = m.find_pairings
puts p.sort.to_s # gives an array of arrays, the first value in each sub-array is the original position, the second is the new position