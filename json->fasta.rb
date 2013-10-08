require "rubygems"
require "json"

def fasta_array (json)
	frags = JSON.parse(File.open(json).read) #open the json containing an array of frags, each frag is an array of nucleotides
	frag_strings = [] #sequences
	id_and_length = [] 
	frag_ids = [] #create an id for each of the frags
	x = 0
	frags.each do |i|
		frag_strings << i.join.to_s #add these to a new array of the sequences in string format
		frag_ids << ('>frag' + (x+1).to_s)
		id_and_length << (('>frag' + (x+1).to_s) + "  Length = " + i.length.to_s)
		x+=1
	end
	fastaformat_array = id_and_length.zip(frag_strings) #create the array, each element of which goes onto a new line in fasta
	return fastaformat_array, frags
end

def vcf_array (frags)
	snp_pos = JSON.parse(File.open('snp_pos.json').read)#open the json containing the position data for the snps on each fragment - array of arrays
	ids = []
	x = 0
	frags.each do |i| #getting the ids - removing the '>'
		ids << ('frag' + (x+1).to_s)
		x+=1
	end
	chrom = []
	alt = []
	q = 0
	snp_pos.each do |h|
		if h.to_s != '[]' #all of the fragments that contain at least one snp, TESTED = 100, so correct
			h.length.times do #*no. of snps
				chrom << ids[q]
			end
		end
		h.each do |i|
			alt << frags[q][i].capitalize #what nucleotide is at these positions?  VCF requires capital nucleotides
		end
		q += 1
	end
	vcf_format = ['##fileformat=VCFv4.1', '##source=Fake', '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'] 
	u = 0
	snp_pos.flatten.each do |i|
		line = chrom[u] + '	' + i.to_s + '	' + '.' + '	' + '.' + '	' + alt[u] + '	' + '100' + '	' + 'PASS' + '	' + '.'
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

fasta_and_frags = fasta_array('frags.json')
fastaformat_array = fasta_and_frags[0]
#fastaformat_array_shuf = fastaformat_array.shuffle #shuffle it for shits and giggles #YOU CAN USE THE SHUFFLED VERSION LATER ON
frags = fasta_and_frags[1]
vcf = vcf_array(frags)
write_fasta(fastaformat_array, 'frags.fasta')
write_vcf(vcf, 'snps.vcf')



