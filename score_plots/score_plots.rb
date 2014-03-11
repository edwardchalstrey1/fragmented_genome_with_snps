#encoding: utf-8
require '~/fragmented_genome_with_snps/lib/write_it'
require '~/fragmented_genome_with_snps/lib/rearrangement_score'
require 'rinruby'

def original_order
	original_order = []
	IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/correct_permutation.txt") { |line| original_order << line }
	original_order[1..-1]
end

def get_perms(gen, start, inc) # number of generations to choose from
	all_perms = []
	x = start
	gen.times do
		pop = []
		Dir.entries("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{x}").each do |ptxt|
			if ptxt.include? '.txt'
				perm = []
				IO.foreach("#{Dir.home}/fragmented_genome_with_snps/arabidopsis_datasets/#{ARGV[0]}/#{ARGV[1]}/Gen#{x}/#{ptxt}") { |line| perm << line }
				pop << perm
			end
		end
		all_perms << pop
		x+=inc # TODO change back to 1
	end
	return all_perms
end

def get_metrics(all_perms)
	orig = original_order
	all_metrics = []
	all_perms.each do |pop|
		pop_met = []
		pop.each do |perm|
			metrics = {}
			metrics[:fit] = (1 - (perm[0].gsub(/\n/, "")).to_f)
			metrics[:dev] = RearrangementScore::dev_dist(orig, perm[1..-1])
			metrics[:sq] = RearrangementScore::sq_dev_dist(orig, perm[1..-1])
			metrics[:ham] = RearrangementScore::gen_ham_dist(orig, perm[1..-1])
			metrics[:mod] = RearrangementScore::mod_ham_dist(orig, perm[1..-1])
			metrics[:r] = RearrangementScore::r_dist(orig, perm[1..-1])
			metrics[:lcs] = RearrangementScore::lcs(orig, perm[1..-1])
			metrics[:kt] = RearrangementScore::kendalls_tau(orig, perm[1..-1])
			pop_met << metrics
		end
		all_metrics << pop_met
	end
	return all_metrics
end

def gg_plots(all_metrics, start, inc)
	myr = RinRuby.new(echo = false)
	myr.eval 'source("~/fragmented_genome_with_snps/score_plots/score_plots.R")'
	x, y, se = [], [], []
	gen = start
	all_metrics.each do |pop_met|
		fit, dev, sq, ham, mod, r, lcs, kt = [],[],[],[],[],[],[],[]
		pop_met.each do |metrics|
			fit << metrics[:fit]
			dev << metrics[:dev]
			sq << metrics[:sq]
			ham << metrics[:ham]
			mod << metrics[:mod]
			r << metrics[:r]
			lcs << metrics[:lcs]
			kt << metrics[:kt]
		end
		pop_y = [fit, dev, sq, ham, mod, r, lcs, kt]
		#pop_y = [dev, sq]
		pop_y.each do |metrics|
			y << metrics.inject(:+) / metrics.length.to_f
			x << gen
			myr.assign 'metrics', metrics
			sem = myr.pull 'st_err(metrics)'
			se << sem
		end
		gen+=inc # TODO change back to 1
	end
	group = ['rev_fit', 'dev', 'sq', 'ham', 'mod', 'rev_r', 'lcs', 'kt'] * all_metrics.length
	#group = ['fit', 'dev'] * (all_metrics.length / 2)
	myr.assign 'x', x
	myr.assign 'y', y
	myr.assign 'se', se
	myr.assign 'group', group
	myr.assign 'dataset_run', "#{ARGV[0]}/#{ARGV[1]}"
	myr.eval 'plot_it(x, y, se, group, dataset_run)'
	myr.quit
end

all_perms = get_perms(3, 0, 125) # [][][0] is the fitness of the permutation [][][1..-1] is the permutation (frag ids)
all_metrics = get_metrics(all_perms)
puts all_metrics[0][0]

gg_plots(all_metrics, 0, 125)

