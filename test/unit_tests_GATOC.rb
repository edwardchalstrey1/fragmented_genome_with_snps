#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/reform_ratio'
require '~/fragmented_genome_with_snps/lib/GATOC'
require '~/fragmented_genome_with_snps/lib/snp_dist'
require '~/fragmented_genome_with_snps/lib/model_genome'
require 'test/unit'

class TestGATOC < Test::Unit::TestCase

	TEST_ARRAY = %w(a b c d e f g h i j k l m n o p q r s t)

	Div = 100.0; Genome_length = 2000.0
	hm, ht = ModelGenome::get_snps('hm <- rnorm(50, 1000, 100)', 'ht <- runif(50, 1, 2000)')
	Fratio, Breaks = SNPdist::fratio(hm, ht, Div, Genome_length)
	Hyp = SNPdist::hyp_snps([Fratio, Breaks], Div, Genome_length)

	Fasta_array = ReformRatio::fasta_array('arabidopsis_datasets/small_dataset2/frags.fasta')
	Snp_data = ReformRatio::get_snp_data('arabidopsis_datasets/small_dataset2/snps.vcf')

	pop = GATOC::initial_population(Fasta_array, 10)
	Selected = GATOC::select(pop, Snp_data, 5, Hyp, Div, Genome_length)

	def test_fitness
		fit, hm, ht, hyp = GATOC::fitness(Fasta_array, Snp_data, 'gen', Hyp, 'loc', 'dat', 'run', Div, Genome_length)
		assert_kind_of(Float, fit)
		assert_in_delta(0.5, fit, 0.5) # 0.5 +- 0.5
		assert_kind_of(Array, hm)
		assert_kind_of(Array, ht)
		assert_kind_of(Array, hyp)
		assert_kind_of(Integer, hm[0])
		assert_kind_of(Integer, ht[0])
		assert_kind_of(Float, hyp[0])
	end

	def test_initial_population
		array = %w(a b)
		pop = GATOC::initial_population(array, 2)
		assert(pop == [%w(a b), %w(a b)] || pop == [%w(b a), %w(a b)] || pop == [%w(a b), %w(b a)] || pop = [%w(b a), %w(b a)])
	end

	def test_select
		assert_kind_of(Integer, Selected[1], 'leftover not int') # leftover
		assert_kind_of(Array, Selected[0], 'permutation and correlation not array') # permutation and correlation
		assert_equal(5, Selected[0].length) # no. of permutations selcted
		assert_kind_of(Float, Selected[0][0][0], 'correlation not float') # correlation value
		assert_in_delta(0.5, Selected[0][0][0], 0.5) # correlation value
		assert_kind_of(Array, Selected[0][0][1], 'permutation not array') # permutation
		assert_kind_of(Bio::FastaFormat, Selected[0][0][1][0], 'Not a Bio::FastaFormat')
		assert_equal(10, Selected[2].length)
	end

	def test_new_population
		new_pop  = GATOC::new_population(Selected[0], 10, 1, 1, 1, 5, Selected[1])
		assert_kind_of(Array, new_pop)
		assert_kind_of(String, new_pop[1])
		assert_kind_of(Array, new_pop[0])
		x = 0
		new_pop[0].each do |permutation, type|
			assert_kind_of(Array, permutation, "permutation at element #{x} of population not array")
			assert_equal(53, permutation.length, "permutation at element #{x} of population not correct number of frags")
			assert_kind_of(String, type, 'type not a string')
			x+=1
		end
		assert_kind_of(Bio::FastaFormat, new_pop[0][0][0][0])
		assert_equal(new_pop.uniq.length, new_pop.length)
	end
end

