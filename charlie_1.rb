require 'rubygems'
require 'bio-samtools'
require 'bio'
require 'rinruby'
require 'charlie'

FRAGS = []
Bio::FastaFormat.open('arabidopsis_datasets/'+ARGV[0].to_s+'/frags_shuffled.fasta').each do |i| #get array of fasta format frags
    FRAGS << i
end

class FragmentArrangement < PermutationGenotype(FRAGS.size)
    def get_snp_data (vcf_file)
        vcfs_chrom = []
        vcfs_pos = []
        vcfs_info = []
        File.open(vcf_file, "r").each do |line| #get array of vcf lines, you can call a method on one line
            next if line =~ /^#/
                v = Bio::DB::Vcf.new(line)
                vcfs_chrom << v.chrom
                vcfs_pos << v.pos
                vcfs_info << v.info # so this will be an array of hashes of strings
            end
            num_snps_frag_hash = Hash.new(0)
            vcfs_chrom.each {|v| num_snps_frag_hash[v] +=1 } #we have the number of snps on each frag, by counting the repeats of each frag in the vcf
            #the frag_id(.chrom) is the key, the number of snps for that frag is the value. putting the number of snps for each frag into hash
            return vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info
    end
    def fasta_array #(fasta_file)
        ids = []
        #fasta = [] #we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
        lengths = []
        #Bio::FastaFormat.open(fasta_file).each do |i| #get array of fasta format frags, ##  WE NEED TO REORDER THE FASTA FRAGS HERE, TO TEST DIFFERENT ARRANGEMENTS
        #        fasta << i
        #end
        #fasta = fasta.shuffle
        #fasta = fasta_child(fasta)
        fasta = FRAGS
        fasta.each do |i|
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
    def total_pos (pos, fasta_lengths) #pos is in the same order as the vcf, fasta lengths is in the fasta order, so this only works if they are in the correct order, 
        totals = []                    #but the totals are what we need for wrongly ordered frags too, they will be incorrect of course, but that is what we need to test
        x = 0                                                   # matching the ids for the info_hash and total_pos will work, because total pos is in the same order as pos, same order as VCF
        pos.each do |frag|
            if x == 0
            	totals << frag.uniq
                x+=1
            else
                tot_frag = []
                lengths = []
                fasta_lengths[0..x-1].each do |p|
                	lengths << p
            	end
            	so_far = lengths.inject(:+) # this needs to be the length of the frags, not the number of snps
            	frag.each do |i|
            	    tot_frag << so_far-1 + i
            	end
            	totals << tot_frag.uniq
            	x+=1
        	end
    	end
        return totals
    end
    def het_hom (actual_pos, vcfs_info) #actual_pos in same order as VCF
        het = []
        hom = []
        x = 0
        actual_pos.flatten.each do |snp|
            if vcfs_info[x] == {"AF"=>"1.0"} # homozygous SNPs have AF= 1.0, ###we can change this to a range for real data###
                hom << snp
            elsif vcfs_info[x] == {"AF"=>"0.5"}
                het << snp
            end
            x+=1
        end
        return het, hom
    end
    def fitness
        snp_data = get_snp_data('arabidopsis_datasets/'+ARGV[0].to_s+'/snps.vcf')
        vcfs_chrom = snp_data[0] #array of vcf frag ids
        vcfs_pos = snp_data[1] #array of all the snp positions (fragments with snps)
        snps_hash = snp_data[2] #hash of each fragment from vcf, and it's number of snps
        vcfs_info = snp_data[3]
        fasta_data = fasta_array #array of fasta format fragments, and entry_ids. 					###	THIS IS IMPORTANT. NEEDS TO BE THE PERMUTATION! ###
        fasta = fasta_data[0]
        fasta_ids = fasta_data[1]
        fasta_lengths = fasta_data[2]
        snps_per_frag = snps_per_fasta_frag(snps_hash, fasta) #array of no. of snps per frag in same order as fasta
        pos = get_positions(fasta, vcfs_chrom, vcfs_pos, snps_per_frag) #get snp positions for each frag in array of arrays
        actual_pos = total_pos(pos, fasta_lengths)
        het_hom_snps = het_hom(actual_pos, vcfs_info)
        het = het_hom_snps[0]
        hom = het_hom_snps[1]
        myr = RinRuby.new(echo=false)
        myr.assign "het_snps", het
        myr.assign "hom_snps", hom
        myr.eval "source('~/fragmented_genome_with_snps/ratio.R')"
        coeff = myr.pull "cor(qqp$x,qqp$y)"
        myr.quit
        return coeff
    end
    use RouletteSelection, EdgeRecombinationCrossover
end
Population.new(FragmentArrangement, 20).evolve_on_console(10)
