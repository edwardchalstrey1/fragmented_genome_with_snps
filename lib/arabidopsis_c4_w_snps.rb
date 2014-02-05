#encoding: utf-8
require 'rubygems'
require 'rinruby'
require 'json'
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
	myr.eval "hm <- rnorm(35, 10000000, 5000000)"
	myr.eval "ht1a <- rnorm(1500, 5000000, 1000000)"
	myr.eval "ht1 <- ht1a[which(ht1a < 7.5e+06)]"  #non-recombinant region = 7.5m-12.5m, don't want heterozygous SNPs here
	myr.eval "ht2a <- rnorm(1500, 15000000, 1000000)"
	myr.eval "ht2 <- ht2a[which(ht2a > 1.25e+07)]"
	myr.eval "ht <- c(ht1, ht2)"
	hm = myr.pull "hm"
	unless hm.include?(10000000)
		hm = [hm, 10000000].flatten # adding in the causative SNP, if there isn't one already at the position
	end
	ht = myr.pull "ht"
	hm_int, ht_int = []
	hm.each{|i| hm_int << i.abs.to_i}
	ht.each{|i| ht_int << i.abs.to_i}
	snp_pos = [ht_int,hm_int].flatten
	myr.quit
	puts snp_pos.uniq.length.to_s+" of "+snp_pos.length.to_s+" SNPs are unique"
	puts "There are "+hm.length.to_s+" homozygous SNPs"
	puts "There are "+ht.length.to_s+" heterozygous SNPs"
	return snp_pos, hm_int, ht_int
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
		case seq[i-1]
		when 'A' then seq[i-1] = 'T'
		when 'T' then seq[i-1] = 'A'
		when 'C' then seq[i-1] = 'G'
		when 'G' then seq[i-1] = 'C'
		when 'N' then seq[i-1] = 'R'
		else puts "Error: base not A,T,C,G or N\n Base was #{seq[i-1]}"
		end
	end
	return seq
end

def pos_each_frag (snp_pos, frags) # get the positions for snps on individual frags
	snp_pos.sort! # this is needed as the ht/hm SNPs need to be ordered together
	p_ranges = []
	frags.each do |i|
		p_ranges << i.length + p_ranges[-1].to_i # adding the fragment lengths to get the upper bounds of ranges of positions on the original seq.
	end
	first_pos = [] # then, to work out the first position of each fragment
	first_pos << 0
	first_pos << p_ranges[0..-2]
	first_pos = first_pos.flatten
	p = 0
	t = 0
	all_frags_pos = [] # the positions of snps on each fragment, array of arrays
	snp_pos_all = [] # the actual positions in the genome for each fragment
	p_ranges.each do |jj| # for each of the upper ranges (j) that the positions could be in
		each_frag_pos = []
		actual_pos = []
		while snp_pos[t].to_i < jj && snp_pos[t].nil? == false do # make the loop quit before trying to access index of snp_pos that doesn't exist
			each_frag_pos << snp_pos[t].to_i 	# add all of the positions < jj to a new array for each frag
			actual_pos << snp_pos[t].to_i
			t += 1
		end
		snp_pos_all << actual_pos 
		z = 0
		y = first_pos[p].to_i
		each_frag_pos.each do |e|
      each_frag_pos[z] = e - y # taking away the value at the first position of each fragment to get the true/actual positions
      z += 1
    end
		p += 1
		all_frags_pos << each_frag_pos # add the positions for the frag to an array of the each_frag arrays
	end
	return all_frags_pos, snp_pos_all
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

def write_txt (filename, array)
	File.open(filename, "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end

###########
#FASTA/VCF#
###########

def fasta_array (frags)
	frag_ids, id_and_length = [] 
	x = 0
	frag_strings = []
	frags.each do |i|
		id = ('>frag' + (x+1).to_s)
		id_and_length << (id + "  Length = " + i.length.to_s)
		id.slice!('>')
		frag_ids << id
		frag_strings << i.join.to_s # joining all the bases in one frag array to form a string for fasta
		x+=1
	end
	fastaformat_array = id_and_length.zip(frag_strings) #create the array, each element of which goes onto a new line in fasta
	return fastaformat_array, frag_ids
end

def vcf_array (frags, pos_on_frags, snp_pos_all, hm, ht)
	x = 0
	ids = []
	frags.each do |i| #getting the ids 
		ids << ('frag' + (x+1).to_s)
		x+=1
	end
	chrom, alt = []
	q = 0
	pos_on_frags.each do |h|
		if h.length != 0 #all of the fragments that contain at least one snp
			h.length.times do #*no. of snps
				chrom << ids[q]
			end
		end
		h.each do |i|
			alt << frags[q][i] #what nucleotide is at these positions?
		end
		q += 1
	end
	ref = []
	it = 0
	alt.each do |base| # TODO change to CASE!!! NOT SURE IF THIS IS RIGHT, WONT IT NOW SKIP THE CASE IF BASE IS nil? I WANT IT TO DO BOTH
		base == nil ? alt[it] ='T': case base
		when 'A' then ref << 'T'
		when 'T' then ref << 'A'
		when 'C' then ref << 'G'
		when 'G' then ref << 'C'
		when 'R' then ref << 'N'
		else ref << 'A'
		end
		it+=1
	end
	vcf_format = ['##fileformat=VCFv4.1', '##source=Fake', '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'] 
	u = 0
	pos_on_frags.flatten.each do |i| # trying to do: if we increment the pos_on_frags and the snp_pos_all together, we should be able to tell whether each SNP from on_frags is het/homo
		if hm.include?(snp_pos_all.flatten[u])
			x = "1.0"
		elsif ht.include?(snp_pos_all.flatten[u])
			x = "0.5"
		else
			x = "WRONG"
		end
		line = chrom[u] + '	' + i.to_s + '	.	' +	ref[u] + '	' + alt[u] + '	100	PASS	AF='+x
		vcf_format << line
		u += 1
	end
	return vcf_format, chrom.uniq
end

def write_data (array, file, dataset)
	File.open("arabidopsis_datasets/"+dataset.to_s+"/"+file, "w+") do |f|
		array.each { |i| f.puts(i) } #write the fasta/vcf
	end
end

def json_array (array, file, dataset)
	File.open("arabidopsis_datasets/"+dataset.to_s+"/"+file, "w+") do |f|
		f.write(array.to_json)
	end
end

Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/"+ARGV[0].to_s)) # make the directory to put data files into

snpz = normal_dist
snp_pos = snpz[0]
hm = snpz[1]
ht = snpz[2]

puts "Is there a SNP at the causative mutation position? -- "+snp_pos.include?(10000000).to_s

arabidopsis_c4 = fasta2char_array("TAIR10_chr4.fasta")
snp_sequence = snp_seq(arabidopsis_c4, snp_pos)
frags = get_frags(snp_sequence)
puts "Arabidopsis chr4 length: "+arabidopsis_c4.length.to_s+" bases"
puts "Fragmented seq   length: "+frags.join.length.to_s+ " = close enough? You decide."
puts "You have created "+frags.length.to_s+" fragments of sizes 10-20Kb"

pos_on_all = pos_each_frag(snp_pos, frags)
pos_on_frags = pos_on_all[0]
snp_pos_all = pos_on_all[1]

fasta_n_ids = fasta_array(frags)
fastaformat_array = fasta_n_ids[0]
frag_ids = fasta_n_ids[1]

vcf_n_chrom = vcf_array(frags, pos_on_frags, snp_pos_all, hm, ht)
vcf = vcf_n_chrom[0]
chrom = vcf_n_chrom[1]

write_data(fastaformat_array, 'frags.fasta', ARGV[0])
write_data(vcf, 'snps.vcf', ARGV[0])

fastaformat_array_shuf = fastaformat_array.shuffle #shuffle it to show that the order doesn't need to be conserved when working out density later on
write_data(fastaformat_array_shuf, 'frags_shuffled.fasta', ARGV[0])

json_array(frag_ids, 'frag_ids_original_order.json', ARGV[0]) #can't remember why this is needed exactly

################################################################ Everything below for skew scatter graph
def skew_prereq (frags, chrom, pos_on_frags)
	Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/"+ARGV[0].to_s+"/skew_scatter"))
	write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/skew_scatter/snp_pos.txt', pos_on_frags)
	lengths = []
	frags.each do |frag| 
		lengths << frag.length
	end
	frags_w_snps = []
	chrom.each do |id|
		id.slice!("frag")
		frags_w_snps << id.to_i
	end
	lengths_fws = []
	frags_w_snps.each do |id|
		write_txt("arabidopsis_datasets/"+ARGV[0].to_s+"/skew_scatter/snps"+id.to_s+".txt", pos_on_frags[id-1]) #need positions and lengths of each fragment for super skew scatter
		lengths_fws << lengths[id-1]
	end
	write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/skew_scatter/ex_fasta_lengths.txt', lengths_fws)
	write_txt('arabidopsis_datasets/'+ARGV[0].to_s+'/skew_scatter/ex_ids_w_snps.txt', frags_w_snps)
end
skew_prereq(frags, chrom, pos_on_frags)
