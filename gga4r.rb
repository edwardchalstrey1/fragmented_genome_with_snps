require "rubygems"
require "gga4r"

class StringPopulation < Array # use the term self to act on the array/StrPop object in question, i.e. the array of fragments
	def fitness
	
	end

	def recombine(c2)
	cross_point = (rand * c2.size).to_i
	c1_a, c1_b = self.separate(cross_point)
	c2_a, c2_b = c2.separate(cross_point)
	StringPopulation.new(c1_a + c2_b)
	end
end

def create_population_with_fit_all_1s(num = 10, frags_array)
	population = []
	num.times do
		chromosome = StringPopulation.new(frags_array) #here we enter the shuffled/random order fragment array
		population << chromosome
	end
	population
end

ga = GeneticAlgorithm.new(create_population_with_fit_all_1s)

5.times { |i|  ga.evolve }
#p ga.best_fit[0]