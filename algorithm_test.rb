require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'
require_relative 'lib/reform_ratio'

myr = RinRuby.new(echo=false)
myr.eval "source('~/fragmented_genome_with_snps/lib/comparable_ratio.R')"
RATIO = myr.pull "comparable_ratio(1)"
myr.quit

vcf_file = "arabidopsis_datasets/#{ARGV[0]}/snps.vcf"
fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags_shuffled.fasta"

ordered_fasta_file = "arabidopsis_datasets/#{ARGV[0]}/frags.fasta"
ordered_fasta = ReformRatio::fasta_array(ordered_fasta_file)

ReformRatio::evolve(fasta_file, vcf_file, 10, 20, 10, 2, 1, 1, ordered_fasta, "no figures")
