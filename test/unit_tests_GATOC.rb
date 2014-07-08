#encoding: utf-8
require_relative '../lib/reform_ratio'
require_relative '../lib/GATOC'
require_relative '../lib/snp_dist'
require_relative '../lib/model_genome'
require 'test/unit'

class TestGATOC < Test::Unit::TestCase

	def setup
		@fasta_array = ReformRatio::fasta_array('arabidopsis_datasets/small_dataset2/frags.fasta')
		@snp_data = ReformRatio::get_snp_data('arabidopsis_datasets/small_dataset2/snps.vcf')
		@genome_length = ReformRatio::genome_length('arabidopsis_datasets/small_dataset2/frags.fasta')
		@pop = GATOC::initial_population(@fasta_array, 10)
		@selected = GATOC::select(@pop, @snp_data, 5, @genome_length)
	end

	def test_fitness
		fit, hm, ht = GATOC::fitness(@fasta_array, @snp_data, @genome_length)
		assert_kind_of(Float, fit)
		assert_kind_of(Array, hm)
		assert_kind_of(Array, ht)
		assert_kind_of(Integer, hm[0])
		assert_kind_of(Integer, ht[0])
	end

	def test_initial_population
		assert_equal(10, @pop.length)
		assert_kind_of(Array, @pop[0])
		assert_kind_of(Bio::FastaFormat, @pop[0][0][0])
		@pop.each do |perm, type|
			assert_equal(53, perm.length)
		end
	end

	def test_select
		assert_kind_of(Integer, @selected[1], 'leftover not int') # leftover
		assert_kind_of(Array, @selected[0], 'permutation and correlation not array') # permutation and correlation
		assert_equal(5, @selected[0].length) # no. of permutations selcted
		assert_kind_of(Float, @selected[0][0][0], 'correlation not float') # correlation value
		assert_kind_of(Array, @selected[0][0][1], 'permutation not array') # permutation
		assert_kind_of(Bio::FastaFormat, @selected[0][0][1][0], 'Not a Bio::FastaFormat')
		assert_equal(10, @selected[2].length)
	end

	def test_new_population
		new_pop = GATOC::new_population(@selected[0], 10, 4, 2, 2, 2, 5, @selected[1])
		assert_kind_of(Array, new_pop)
		assert_kind_of(Bio::FastaFormat, new_pop[0][0][0])
		assert_equal(new_pop.uniq.length, new_pop.length)
		assert_equal(10, new_pop.length)
		x = 0
		new_pop.each do |permutation, type|
			assert_kind_of(Array, permutation, "permutation at element #{x} of population not array")
			assert_equal(53, permutation.length, "permutation at element #{x} of population not correct number of frags")
			assert_kind_of(String, type, 'type not a string')
			x+=1
		end
	end

	def test_quit
		auc = GATOC.quit([0.9,0.9,0.92,0.95,0.955])
		assert_kind_of(Float, auc)
		assert_equal(3.6975, ('%.4f' % auc).to_f)
	end
end

