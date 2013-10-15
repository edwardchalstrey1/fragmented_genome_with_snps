require 'rubygems'
require 'bio-samtools'
require 'bio'
require "json"

def odds (array) #create array from odd indexes of input array
	return array.values_at(*array.each_index.select {|i| i.odd?})
end
def evens (array)
	return array.values_at(*array.each_index.select {|i| i.even?})
end
def get_snp_data (vcf_file)
	vcfs_chrom = []
	vcfs_pos = []
	File.open(vcf_file, "r").each do |line| #get array of vcf lines, you can call a method on one line
		next if line =~ /^#/
		v = Bio::DB::Vcf.new(line)
		vcfs_chrom << v.chrom
		vcfs_pos << v.pos
	end
	num_snps_frag_hash = Hash.new(0)
	vcfs_chrom.each {|v| num_snps_frag_hash[v] +=1 } #we have the number of snps on each frag, by counting the repeats of each frag in the vcf
	#the frag_id(.chrom) is the key, the number of snps for that frag is the value. putting the number of snps for each frag into hash
	return vcfs_chrom, vcfs_pos, num_snps_frag_hash
end
def fasta_array (fasta_file)
	ids = []
	fasta = [] #we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
	lengths = []
	Bio::FastaFormat.open(fasta_file).each do |i| #get array of fasta format frags
		fasta << i
		ids << i.entry_id
		lengths << i.length
	end
	return fasta, ids, lengths
end
def snps_per_fasta_frag (snps_per_vcf_frag_hash, fasta_array)
	snps_per_frag_fasta_order = [] #use the id to identify the number of snps for that frag using the keys of snps_hash
	fasta_array.each do |frag|
		snps_per_frag_fasta_order << snps_per_vcf_frag_hash[frag.entry_id] #gives 0 for keys that don't exist = good, because the frags with 0 density would otherwise be missing
	end
	#now we have an array with the number of snps per frag in the same order as the fasta array, which we can get lengths from to calculate density
	return snps_per_frag_fasta_order
end
def calculate_densities (fasta, snps_per_frag) #argument arrays must be in "same order" - each fasta frag corresponds to the equivalent number of snps at the same index in snps_per_frag
	densities = []
	x = 0
	fasta.each do |frag|
		if snps_per_frag[x] == 0 #snps_per_frag and fasta in the same order
			densities << 0
			x+=1
		else snps_per_frag[x]
			densities << (snps_per_frag[x].to_f / frag.length.to_f)*1000 #this gives snps/Kb as density units
			x+=1
		end
	end
	return densities
end
def density_order_ids (densities, fasta_ids) #must be in same order
	density_order = (densities.zip(fasta_ids)).sort # array of arrays (density and id values) - these will of course be in the same order
	#the sort, sorts by the first value of each sub array (density)
	frags_by_density = odds(density_order.flatten)
	sorted_densities = evens(density_order.flatten)
	# selecting the odd values from the flattened density order array produces a list of frag ids ranked by density
	# the frags with zero snps are sorted by their id, this doesn't matter, we assume they are in a random unknowable order
	return frags_by_density, sorted_densities
end
def get_positions (fasta, vcfs_chrom, vcfs_pos, snps_per_frag)
	pos = [] #get the snp positions for each frag, in an array of arrays
	n = 0
	fasta.each do |frag|
		x = 0
		each_fr_pos = []
		snps_per_frag[n].times do |j|
			if frag.entry_id == vcfs_chrom[x] #this assumes that frag_id == vcf.chrom then continues to for the number of snps (for that frag)
				each_fr_pos << vcfs_pos[x]
				x+=1
			else
				while frag.entry_id != vcfs_chrom[x]
					x+=1
				end
				each_fr_pos << vcfs_pos[x]
				x+=1
			end
		end
		pos << each_fr_pos #this gives empty arrays for frags with out snps, and a list of the positions of those with
		n+=1
	end
	return pos
end
def associate_fasta_ids_snp_pos (pos, fasta_ids)
	#make hash where key is frag_id and values are positions for each frag, positions array should be converted to string so they stick together
	#otherwise we have problems when flattening the array to make the hash below
	pos_strings = []
	pos.each do |s|
		pos_strings << "#{s.join(",")}" #each position separated by a comma
	end
	return Hash[*fasta_ids.zip(pos_strings).flatten] # fasta ids as keys and the snp position strings as values
end
def write_json (array, json)
	File.open(json, "w") do |f|
		f.write(array.to_json)
	end
end
def two_a_to_h (keys_array, values_array)
	return Hash[*keys_array.zip(values_array).flatten]
end

snp_data = get_snp_data('snps.vcf')
vcfs_chrom = snp_data[0] #array of vcf frag ids
vcfs_pos = snp_data[1] #array of all the snp positions (fragments with snps)
snps_hash = snp_data[2] #hash of each fragment from vcf, and it's number of snps

fasta_data = fasta_array('frags_shuffled.fasta') #array of fasta format fragments, and entry_ids
fasta = fasta_data[0]
fasta_ids = fasta_data[1]
fasta_lengths = fasta_data[2]

snps_per_frag = snps_per_fasta_frag(snps_hash, fasta) #array of no. of snps per frag in same order as fasta

densities = calculate_densities(fasta, snps_per_frag)

sorted_frags_densities = density_order_ids(densities, fasta_ids)
frags_by_density = sorted_frags_densities[0]
write_json(frags_by_density, 'frags_by_density.json')


pos = get_positions(fasta, vcfs_chrom, vcfs_pos, snps_per_frag) #get snp positions for each frag in array of arrays

pos_hash = associate_fasta_ids_snp_pos(pos, fasta_ids) #associate frag ids with a string of the positions it has THIS IS FOR THE LEFT/RIGHT METHOD
write_json(pos_hash, 'pos_hash.json')
write_json(fasta_lengths, 'fasta_lengths.json') #the lengths of the frags in the same order as the fasta file and hence pos_hash

id_density_hash = two_a_to_h(frags_by_density, sorted_frags_densities[1]) #the frags with associated densities sorted
write_json(id_density_hash, 'id_density_hash.json')
