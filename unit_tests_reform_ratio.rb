#encoding: utf-8
require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'
require_relative 'lib/reform_ratio'
require 'test/unit'

class FakeFasta
	attr_accessor :entry_id, :length
	def initialize
		@entry_id = 'FragX'
		@length = 10
	end
end

class TestReform < Test::Unit::TestCase

	TEST_ARRAY = %w(a b c d e f g h i j k l m n o p q r s t)

	def fasta
		frag1, frag2, frag3 = FakeFasta.new, FakeFasta.new, FakeFasta.new
		frag1.entry_id, frag2.entry_id, frag3.entry_id = 'frag1', 'frag2', 'frag3'
		[frag1, frag2, frag3]
	end

	def test_rearrangement_score
		a = %w(a b c)
		c = a.reverse
		assert_equal(0, ReformRatio::rearrangement_score(a,a))
		assert_equal(4, ReformRatio::rearrangement_score(c,a))
	end

	def test_total_pos
		pos = [[1,2,3],[4,5,6],[7,8,9]]
		lengths = [10,10,10]
		expected = [1,2,3,13,14,15,26,27,28]
		assert_equal(expected, ReformRatio::total_pos(pos, lengths))
	end

	def test_prime?
		assert_equal(true, ReformRatio::prime?(7))
		assert_equal(false, ReformRatio::prime?(4))
		assert_equal(true, ReformRatio::prime?(2))
		assert_equal(true, ReformRatio::prime?(971))
		assert_equal(false, ReformRatio::prime?(972))
	end

	def test_division
		a = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
		b = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
		ax = ReformRatio::division(a)
		bx = ReformRatio::division(a)
		assert(ax==2||ax==1||ax==4||ax==5||ax=10||ax==20) 
		assert(bx==18||bx==9||bx=6||bx==3||bx==2||bx==1)
	end

	def test_recombination
		parent1 = TEST_ARRAY #20
		parent2 = parent1.reverse
		child = ReformRatio::recombine(parent1, parent2)

		parent3 = TEST_ARRAY[0..-2] #19
		parent4 = parent3.reverse
		child2 = ReformRatio::recombine(parent3, parent4)
		
		assert(child.uniq == child, 'Child of p1/2 not unique')
		assert(child != parent1, 'Child same as parent1')
		assert(child != parent2, 'Child same as parent2')

		assert(child2.uniq == child2, 'Child of p3/4 not unique')
		assert(child2 != parent3, 'Child same as parent3')
		assert(child2 != parent4, 'Child same as parent4')
	end

	def test_fasta_id_n_lengths
		ids_n_lengths = ReformRatio::fasta_id_n_lengths(fasta)
		ids = ids_n_lengths[0]
		lengths = ids_n_lengths[1]
		assert_equal(['frag1', 'frag2', 'frag3'], ids)
		assert_equal([10,10,10], lengths)
	end

	def test_snps_per_fasta_frag
		h = {'frag2'=>3,'frag1'=>2,'frag3'=>5}
		assert_equal([2,3,5], ReformRatio::snps_per_fasta_frag(h, fasta))
	end

	def test_get_positions
		vcfs_chrom = ['frag3', 'frag3', 'frag3', 'frag1', 'frag1', 'frag2']
		vcfs_pos = [2,3,5,4,7,5] # f3= 2,3,5  f1= 4,7  f2= 5
		vcfs_info = [{'AF'=>'snp1'},{'AF'=>'snp2'},{'AF'=>'snp3'},{'AF'=>'snp4'},{'AF'=>'snp5'},{'AF'=>'snp6'}] # f3= s1,2,3  f1= s4,5  f2= s6
		snps_per_frag = [2,1,3] # in same order as fasta array
		pos_n_info = ReformRatio::get_positions(fasta, vcfs_chrom, vcfs_pos, snps_per_frag, vcfs_info)
		pos = pos_n_info[0]
		info = pos_n_info[1]
		assert_equal([[4,7],[5],[2,3,5]], pos)
		assert_equal([[{'AF'=>'snp4'},{'AF'=>'snp5'}] ,[{'AF'=>'snp6'}], [{'AF'=>'snp1'},{'AF'=>'snp2'},{'AF'=>'snp3'}]], info)
	end

	def test_fasta_array
		fasta_array = ReformRatio::fasta_array('test/dummy.fasta')
		assert_equal('frag1', fasta_array[0].entry_id)
		assert_equal('AAAAAAAA', fasta_array[1].seq)
		assert_equal(8, fasta_array[2].length)
	end

	def test_get_snp_data
		vcfs_chrom = %w(frag1 frag1 frag2 frag3)
		vcfs_pos = [7,8,2,2]
		num_snps_frag_hash = {'frag1'=>2, 'frag2'=>1, 'frag3'=>1}
		vcfs_info = [{'AF'=>'1.0','NS'=>'5'}, {'AF'=>'1.0'}, {'AF'=>'0.5'}, {'AF'=>'1.0'}]
		snp_data = [vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info]
		assert_equal(snp_data, ReformRatio::get_snp_data('test/dummy.vcf'))
	end

	def test_het_hom
		vcfs_info = [{'AF'=>'1.0'}, {'AF'=>'1.0'}, {'AF'=>'0.5'}, {'AF'=>'0.5'}, {'AF'=>'1.0'}]
		actual_pos = [2, 17, 56, 190, 191]
		hom = [2, 17, 191]
		het = [56, 190]
		assert_equal([het,hom], ReformRatio::het_hom(actual_pos, vcfs_info))
	end

	def test_mutate
		mutant = ReformRatio::mutate(TEST_ARRAY)
		assert(mutant.uniq == mutant)
		assert(mutant != TEST_ARRAY, 'Mutant was the same as parent')
		assert_kind_of(Array, mutant, 'Mutant not an array!')	
	end

	def test_mini_mutate
		mini_mutant = ReformRatio::mini_mutate(TEST_ARRAY)
		assert(mini_mutant.uniq == mini_mutant)
		assert(mini_mutant != TEST_ARRAY, "mini_mutant same as non-mutant")
		assert_kind_of(Array, mini_mutant)
	end

	def test_fitness
		fasta_array = ReformRatio::fasta_array('arabidopsis_datasets/ratio_dataset3/frags.fasta')
		snp_data = ReformRatio::get_snp_data('arabidopsis_datasets/ratio_dataset3/snps.vcf')
		fit = ReformRatio::fitness(fasta_array, snp_data, "diff")
		assert_kind_of(Float, fit)
		assert_in_delta(0.5, fit, 0.5)
	end

	def test_initial_population
		array = %w(a b)
		pop = ReformRatio::initial_population(array, 2)
		assert(pop == [%w(a b), %w(a b)] || pop == [%w(b a), %w(a b)] || pop == [%w(a b), %w(b a)] || pop = [%w(b a), %w(b a)])
	end

	def test_select
		fasta_array = ReformRatio::fasta_array('arabidopsis_datasets/ratio_dataset4/frags_shuffled.fasta')
		snp_data = ReformRatio::get_snp_data('arabidopsis_datasets/ratio_dataset4/snps.vcf')
		pop = ReformRatio::initial_population(fasta_array, 20)
		selected = ReformRatio::select(pop, snp_data, 10)
		assert_kind_of(Integer, selected[1], 'leftover not int') # leftover
		assert_kind_of(Array, selected[0], 'permutation and correlation not array') # permutation and correlation
		assert_equal(10, selected[0].length) # no. of permutations selcted
		assert_kind_of(Float, selected[0][0][0], 'correlation not float') # correlation value
		assert_in_delta(0.5, selected[0][0][0], 0.5) # correlation value
		assert_kind_of(Array, selected[0][0][1], 'permutation not array') # permutation
		assert_kind_of(Bio::FastaFormat, selected[0][0][1][0], 'Not a Bio::FastaFormat')
	end
end

