#encoding: utf-8
class Permutation

	def initialize(name, fitness_score, type, fasta_ids)
		@name = name
		@fitness_score = fitness_score
		@type = type
		@fasta_ids = fasta_ids
	end

	def fitness
		@fitness_score
	end

	def name
		@name
	end

	def type
		@type
	end

	def fasta_ids
		@fasta_ids
	end

end
