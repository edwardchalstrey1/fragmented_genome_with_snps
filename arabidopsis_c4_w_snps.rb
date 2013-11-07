require "rubygems"
require "rinruby"
require "json"
require 'bio-samtools'
require 'bio'

def normal_dist
	myr = RinRuby.new(echo = false)
	myr.eval "x <- rnorm(70000, 10000000, 2000000)"
	snp_pos = myr.pull "x"
	return snp_pos
end
def get_frags (seq)
	seq = seq
	frags = []
	rt = 0
	while rt < seq.length
		frag_length = rand(20000) + 5000
		frag = seq[rt..(rt+frag_length)]
		rt = rt+frag_length
		frags << frag
	end
	return frags
end

snp_pos = normal_dist
puts snp_pos.uniq.length.to_s+" of 70K SNPs are unique"

arabidopsis_c4=[]
arabidopsis_c4<<Bio::FastaFormat.open("TAIR10_chr4.fasta")

#frags = get_frags(arabidopsis_c4)

puts arabidopsis_c4[0].seq #GET THIS TO WORK!!!!
