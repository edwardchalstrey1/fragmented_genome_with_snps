require '~/fragmented_genome_with_snps/lib/write_it'
require 'rinruby'
require 'pp'

# ARGV[0] = dataset, 1 = run, 2 = number of generations to plot

class ReproductionMethods

	def self.get_generations(location, dataset, run, gen_num)
		type_n_fits = [] 
		Dir.chdir(File.join(Dir.home, location)) do
			x = 0
			gen_num.to_i.times do
				type_n_fits << WriteIt::file_to_array("#{dataset}/#{run}/Gen#{x}_lists/gen_#{x}_types.txt")
				x+=1
			end
		end
		x = 0
		type_n_fits.length.times do
			type_n_fits[x] = type_n_fits[x].each_slice(2).to_a
			x+=1
		end
		recomb, mut, minmut, rand, save, x_gen, st_errs, gen = [], [], [], [], [], [], [], 0

		type_n_fits.each do |t_gen| # each generation/population
			re, mu, mi, ra, sa = [], [], [], [], []
			t_gen.each do |type| # permutation
				case type[0]
					when 'recombined' then re << type[1].to_f
					when 'mutant' then mu << type[1].to_f
					when 'mini_mutant' then mi << type[1].to_f
					when 'random' || 'exta_rand' then ra << type[1].to_f
					when 'saved' then sa << type[1].to_f
				end
			end
			if gen > 0
				recomb << re.inject(:+) / re.length.to_f
				mut << mu.inject(:+) / mu.length.to_f
				minmut << mi.inject(:+) / mi.length.to_f
				save << sa.inject(:+) / sa.length.to_f
			end
			rand << ra.inject(:+) / ra.length.to_f

			if gen > 0
				5.times{x_gen << gen}
			else
				x_gen << gen
			end
			puts gen
			gen+=1
			[re, mu, mi, ra, sa].each do |type| # all the fitness' for the permuations of a specific type in one generation (calculate standard error)
				if type.length > 0
					myr = RinRuby.new(echo = false)
					myr.eval 'source("~/fragmented_genome_with_snps/score_plots/reproduction_methods.R")'
					myr.assign 'type', type # array of floats (fitness scores)
					if type.length > 1
						se = myr.pull 'st_err(type)'
					elsif type.length == 1
						se = 0
					end
					myr.quit
					st_errs << se
				end
			end
		end
		y_fits = [recomb, mut, minmut, rand, save].flatten
		gen_0_groups = ['rand']		
		groups = ['recomb', 'mut', 'minmut', 'rand', 'save'] * (gen_num.to_i - 1)
		groups = [gen_0_groups, groups].flatten
		return type_n_fits, y_fits, x_gen, groups, st_errs
	end
end

type_n_fits, y_fits, x_gen, groups, st_errs = ReproductionMethods::get_generations('fragmented_genome_with_snps/arabidopsis_datasets', ARGV[0], ARGV[1], ARGV[2])

puts type_n_fits.length # should be 22
puts
puts y_fits.length
puts
puts x_gen.length
puts
puts groups.length
puts
puts st_errs.length

