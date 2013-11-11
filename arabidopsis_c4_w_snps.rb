require "rubygems"
require "rinruby"
require "json"
require 'bio-samtools'
require 'bio'

def fasta2char_array (fasta)
	fasta_array=[]
	Bio::FastaFormat.open(fasta).each do |i|
		fasta_array << i.seq
	end
	return fasta_array[0].split(//)
end
def normal_dist
	myr = RinRuby.new(echo = false)
	myr.eval "x <- rnorm(70000, 10000000, 2000000)" #distibution about the mean, midpoint of the sequence 
	snp_pos = myr.pull "x"
	snp_pos << 10000000 # there must be a SNP at the causative mutation location
	return snp_pos
end
def get_frags (seq)
	frags = []
	rt = 0
	while rt < seq.length
		frag_length = rand(10000) + 10000 # fragments of 10kb - 20kb
		frag = seq[rt..(rt+frag_length)]
		rt = rt+frag_length
		frags << frag
	end
	return frags
end
def snp_seq (seq, snp_pos)
	snp_pos.each do |i|
		if seq[i] == 'A'
			seq[i] = 'T'
		elsif seq[i] == 'T'
			seq[i] = 'A'
		elsif seq[i] == 'C'
			seq[i] = 'G'
		elsif seq[i] == 'G'
			seq[i] = 'C'
		elsif seq[i] == 'N'
			seq[i] = 'R'
		end
	end
	return seq
end
def pos_each_frag (snp_pos, frags) #get the positions for snps on individual frags
	sorted_pos = snp_pos.sort #this is needed as we have added the mutant SNP on the end
	p_ranges = []
	frags.each do |i|
		p_ranges << i.length + p_ranges[-1].to_i #adding the fragment lengths to get the upper bounds of ranges of positions on the original seq.
	end
	first_pos = [] #then, to work out the first position of each fragment
	first_pos << 0
	first_pos << p_ranges[0..-2]
	first_pos = first_pos.flatten
	p = 0
	t = 0
	all_frags_pos = [] # the positions of snps on each fragment, array of arrays
	p_ranges.each do |jj| #for each of the upper ranges (j) that the positions could be in
		each_frag_pos = []
		while sorted_pos[t].to_i < jj && !sorted_pos[t].nil? do #make the loop quit before trying to access index of snp_pos that doesn't exist
			each_frag_pos << sorted_pos[t].to_i 	# add all of the positions < jj to a new array for each frag
			t += 1
		end
		z = 0
		y = first_pos[p].to_i
		each_frag_pos.each do |e|   
			each_frag_pos[z] = e - y #taking away the value at the first position of each fragment to get the true positions
			z += 1
		end
		p += 1
		all_frags_pos << each_frag_pos # add the positions for the frag to an array of the each_frag arrays
	end
	return all_frags_pos
end
def storage_json (frags, pos, dataset)
	frags_with_positions = []
	x = 0
	frags.each do |i|
		frag_pos_hash = {}
		frag_pos_hash[:frag] = i.join.to_s
		frag_pos_hash[:pos] = pos[x]
		frags_with_positions << frag_pos_hash 
		x+=1
	end
	File.open("arabidopsis_datasets/"+dataset.to_s+"/"+dataset.to_s+"_fwp.json", "w") do |f|
		f.write(frags_with_positions.to_json)
	end
end

snp_pos = normal_dist
puts snp_pos.uniq.length.to_s+" of 70,001 SNPs are unique"
arabidopsis_c4 = fasta2char_array("TAIR10_chr4.fasta")
snp_sequence = snp_seq(arabidopsis_c4, snp_pos)
frags = get_frags(snp_sequence)
puts "Arabidopsis chr4 length: "+arabidopsis_c4.length.to_s+" bases"
puts "Fragmented seq   length: "+frags.join.length.to_s+ " = close enough? You decide."
puts "You have created "+frags.length.to_s+" fragments of sizes 10-20Kb"
pos_on_frags = pos_each_frag(snp_pos, frags)
Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/"+ARGV[0].to_s))
storage_json(frags, pos_on_frags, ARGV[0])

