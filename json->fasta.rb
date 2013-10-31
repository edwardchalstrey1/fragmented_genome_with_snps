require "rubygems"
require "json"

def extract_json (json)
	snp_pos = []
	frags = []
	JSON.parse(File.open(json).read).each do |hash|
		frags << hash["frag"]
		snp_pos << hash["pos"]
	end
	return frags, snp_pos
end
def fasta_array (frags)
	frag_ids = []
	id_and_length = [] 
	x = 0
	frags.each do |i|
		id = ('>frag' + (x+1).to_s)
		id_and_length << (id + "  Length = " + i.length.to_s)
		id.slice!('>')
		frag_ids << id
		x+=1
	end
	fastaformat_array = id_and_length.zip(frags) #create the array, each element of which goes onto a new line in fasta
	return fastaformat_array, frag_ids
end
def vcf_array (frags, snp_pos)
	x = 0
	ids = []
	frags.each do |i| #getting the ids 
		ids << ('frag' + (x+1).to_s)
		x+=1
	end
	chrom = []
	alt = []
	q = 0
	snp_pos.each do |h|
		if h.to_s != '[]' #all of the fragments that contain at least one snp
			h.length.times do #*no. of snps
				chrom << ids[q]
			end
		end
		h.each do |i|
			alt << frags[q][i].capitalize #what nucleotide is at these positions?  VCF requires capital nucleotides
		end
		q += 1
	end
	ref = []
	alt.each do |base|
		if base == 'A'
			ref << 'T'
		else
			ref << 'A'
		end
	end
	vcf_format = ['##fileformat=VCFv4.1', '##source=Fake', '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'] 
	u = 0
	snp_pos.flatten.each do |i|
		line = chrom[u] + '	' + i.to_s + '	.	' +	ref[u] + '	' + alt[u] + '	100	PASS	.'
		vcf_format << line
		u += 1
	end
	return vcf_format
end
def write_fasta (array, file)
	File.open(file, "w+") do |f|
		array.each { |i| f.puts(i) } #write the fasta
	end
end
def write_vcf (array, file)
	File.open(file, "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end
def write_json (array, json)
	File.open(json, "w") do |f|
		f.write(array.to_json)
	end
end
def write_txt (filename, array)
	File.open(filename, "w+") do |f|
		array.each { |i| f.puts(i) }
	end
end

data = extract_json('frags_with_positions.json')
frags = data[0]
snp_pos = data[1]

fasta_n_ids = fasta_array(frags)
fastaformat_array = fasta_n_ids[0]
frag_ids = fasta_n_ids[1]

vcf = vcf_array(frags, snp_pos)

write_fasta(fastaformat_array, 'frags.fasta')
write_vcf(vcf, 'snps.vcf')

fastaformat_array_shuf = fastaformat_array.shuffle #shuffle it to show that the order doesn't need to be conserved when working out density later on
write_fasta(fastaformat_array_shuf, 'frags_shuffled.fasta')

write_json(frag_ids, 'frag_ids_original_order.json')

################### method for getting 0's and 1's (and lengths)
def snp_y_n (frags, snp_pos)
	lengths = []
	each_nuc_snp = []
	x = 0
	frags.each do |frag| # for each base (of each frag) we want to add a zero or a 1, a 1 if this base is a snp
		lengths << frag.length
		nuc_snp = []
		if snp_pos[x] != nil
			y = 0
			until snp_pos[x][y] == nil
				(snp_pos[x][y]-1-nuc_snp.length).times {nuc_snp << 0}
				nuc_snp << 1
				y+=1
			end
			(frag.length-nuc_snp.length).times {nuc_snp << 0}
		else
			frag.length.times {nuc_snp << 0} # for fragments with no snps, each nucleotide has the zero value
		end
		each_nuc_snp << nuc_snp
		x+=1
	end
	return lengths, each_nuc_snp
end

snp_n_lengths = snp_y_n(frags, snp_pos)
lengths = snp_n_lengths[0]
each_nuc_snp = snp_n_lengths[1]

Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/skew_scatter"))
x=1
each_nuc_snp.each do |frag|
	write_txt('skew_scatter/snps'+x.to_s+'.txt', frag) #need positions and lengths of each fragment for super skew scatter
	x+=1
end
write_txt('skew_scatter/ex_fasta_lengths.txt', lengths)