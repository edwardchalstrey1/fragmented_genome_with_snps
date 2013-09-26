#!/usr/bin/ruby

require "rubygems"
require "json"

json = File.open('frags.json').read #open the json containing an array of frags, each frag is an array of nucleotides
frags = JSON.parse(json)
frag_strings = [] #sequences
frags.each do |i|
	j = i.join.to_s #create a string of each sequence
	frag_strings << j #add these to a new array of the sequences in string format
end
frag_ids = [] #create an id for each of the frags
x = 0
frag_strings.each do |f| 
	frag_ids << ('>frag' + (x += 1).to_s)
end
frag_lengths = [] #lengths
frags.each do |l|
	frag_lengths << l.length
end
id_and_length = [] #id and length wanted for description lines of fasta
q = 0
frag_ids.each do |h|
	one = h + "  Length = " + frag_lengths[q].to_s
	id_and_length << one
	q += 1
end
fastaformat_array = id_and_length.zip(frag_strings) #create the array, each element of which goes onto a new line in fasta
#fastaformat_array_shuf = fastaformat_array.shuffle #shiffle it for shits and giggles
File.open("frags.fasta", "w+") do |f|
	fastaformat_array.each { |i| f.puts(i) } #write the fasta
end

json2 = File.open('snp_pos.json').read #open the json containing the position data for the snps on each fragment - array of arrays
snp_pos = JSON.parse(json2)
ids = []
frag_ids.each do |i| #getting the ids - removing the '>'
	i.slice! '>'
	ids << i
end
####VCF fields
chrom = []
q = 0
snp_pos.each do |h|
	if h.to_s != '[]' #all of the fragments that contain at least one snp, TESTED = 100, so correct
		h.length.times do
			chrom << ids[q]#*no. of snps?
		end
	end
	q += 1
end
id = []
snp_pos.each do |i| #id's blank, not real sequences
	c = "."
	id << c
end
ref = []
snp_pos.each do |i| #no reference sequence
	c = "."
	ref << c
end
alt = []
z = 0
snp_pos.each do |j| #the snps (alt)
	one = []
	j.each do |i|
		one << frags[z][i].capitalize #what nucleotide is at these positions?  VCF requires capital nucleotides
	end
	alt << one
	z += 1
end
qual = [] # the quality can just all be 100
snp_pos.each do |i|
	c = 100
	qual << c
end
filter = [] #all pass filter (obviously theres no filter)
snp_pos.each do |i|
	c = "PASS"
	filter << c
end
info = [] #no info
snp_pos.each do |i|
	c = '.'
	info << c
end
puts frags.length
puts chrom.length

vcf_format = ['##fileformat=VCFv4.1', '##source=Fake', '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'] 
#each element needs to be a line, including a value from each of the fields above. DOUBLE CHECK THAT THE FIELDS ARE ALIGNED!
u = 0
snp_pos.flatten.each do |i|
	line = chrom[u] + '	' + i.to_s + '	' + id[u] + '	' + ref[u] + '	' + alt.flatten[u] + '	' + qual[u].to_s + '	' + filter[u] + '	' + info[u]
	vcf_format << line
	u += 1
end
File.open("snps.vcf", "w+") do |f|
	vcf_format.each { |i| f.puts(i) }
end