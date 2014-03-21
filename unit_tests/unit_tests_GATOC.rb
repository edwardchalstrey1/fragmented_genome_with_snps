#encoding: utf-8
require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'
require '~/fragmented_genome_with_snps/lib/reform_ratio'
require '~/fragmented_genome_with_snps/lib/GATOC'
require 'test/unit'

class TestGATOC < Test::Unit::TestCase

	TEST_ARRAY = %w(a b c d e f g h i j k l m n o p q r s t)

	def test_division
		a = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
		b = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
		ax = GATOC::division(a)
		bx = GATOC::division(a)
		assert(ax==2||ax==1||ax==4||ax==5||ax==10) 
		assert(bx==9||bx=6||bx==3||bx==2||bx==1)
	end

	def test_recombination
		parent1 = TEST_ARRAY #20
		parent2 = parent1.reverse
		child = GATOC::recombine(parent1, parent2)

		parent3 = TEST_ARRAY[0..-2] #19
		parent4 = parent3.shuffle
		child2 = GATOC::recombine(parent3, parent4)

		fasta_array = ReformRatio::fasta_array('test/frags_shuffled.fasta').shuffle
		p2_fasta = fasta_array.shuffle
		child3 = GATOC::recombine(fasta_array, p2_fasta)
		
		assert(child.uniq == child, 'Child of p1/2 not unique')
		assert(child != parent1, 'Child same as parent1')
		assert(child != parent2, 'Child same as parent2')
		assert(child.sort == parent1.sort, 'Child not a permutation')

		assert(child2.uniq == child2, 'Child of p3/4 not unique')
		assert(child2 != parent3, 'Child same as parent3')
		assert(child2 != parent4, 'Child same as parent4')
		assert(child2.sort == parent3.sort, 'Child not a permutation')

		assert(child3.uniq == child3, 'Child of fasta not unique')
		assert(child3 != fasta_array, 'Child same as fasta')
		assert(child3 != p2_fasta, 'Child same as fasta2')
	end

	def test_mutate
		mutant = GATOC::mutate(TEST_ARRAY)
		assert(mutant.uniq == mutant)
		assert(mutant != TEST_ARRAY, 'Mutant was the same as parent')
		assert_kind_of(Array, mutant, 'Mutant not an array!')	
	end

	def test_mini_mutate
		mini_mutant = GATOC::mini_mutate(TEST_ARRAY)
		assert(mini_mutant.uniq == mini_mutant)
		assert(mini_mutant != TEST_ARRAY, "mini_mutant same as non-mutant")
		assert_kind_of(Array, mini_mutant)
	end

	def test_fitness
		fasta_array = ReformRatio::fasta_array('test/frags.fasta')
		snp_data = ReformRatio::get_snp_data('test/snps.vcf')
		fit = GATOC::fitness(fasta_array, snp_data, "diff")
		assert_kind_of(Float, fit)
		assert_in_delta(0.5, fit, 0.5)
	end

	def test_average_fitness
		fasta_array = ReformRatio::fasta_array('test/frags.fasta')
		snp_data = ReformRatio::get_snp_data('test/snps.vcf')
		fit = GATOC::average_fitness(fasta_array, snp_data, 5)
		assert_kind_of(Float, fit)
		assert_in_delta(0.5, fit, 0.5)
	end

	def test_initial_population
		array = %w(a b)
		pop = GATOC::initial_population(array, 2)
		assert(pop == [%w(a b), %w(a b)] || pop == [%w(b a), %w(a b)] || pop == [%w(a b), %w(b a)] || pop = [%w(b a), %w(b a)])
	end

	def test_select
		fasta_array = ReformRatio::fasta_array('test/frags_shuffled.fasta')
		snp_data = ReformRatio::get_snp_data('test/snps.vcf')
		pop = GATOC::initial_population(fasta_array, 20)
		selected = GATOC::select(pop, snp_data, 10)
		assert_kind_of(Integer, selected[1], 'leftover not int') # leftover
		assert_kind_of(Array, selected[0], 'permutation and correlation not array') # permutation and correlation
		assert_equal(10, selected[0].length) # no. of permutations selcted
		assert_kind_of(Float, selected[0][0][0], 'correlation not float') # correlation value
		assert_in_delta(0.5, selected[0][0][0], 0.5) # correlation value
		assert_kind_of(Array, selected[0][0][1], 'permutation not array') # permutation
		assert_kind_of(Bio::FastaFormat, selected[0][0][1][0], 'Not a Bio::FastaFormat')
	end

	def test_new_population
		fasta_array = ReformRatio::fasta_array('test/frags_shuffled.fasta')
		snp_data = ReformRatio::get_snp_data('test/snps.vcf')
		pop = GATOC::initial_population(fasta_array, 20)
		selected = GATOC::select(pop, snp_data, 10)
		new_pop  = GATOC::new_population(selected[0], 20, 1, 1, 1, 10, selected[1])
		assert_kind_of(Array, new_pop)
		assert_kind_of(String, new_pop[1])
		assert_kind_of(Array, new_pop[0])
		x = 0
		new_pop[0].each do |permutation|
			assert_kind_of(Array, permutation, "permutation at element #{x} of population not array")
			assert_equal(1239, permutation.length, "permutation at element #{x} of population not correct number of frags")
			x+=1
		end
		assert_kind_of(Bio::FastaFormat, new_pop[0][0][0])
		assert_equal(new_pop.uniq.length, new_pop.length)
	end
end
