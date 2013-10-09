require 'rubygems'
require 'bio-samtools'
require 'bio'
require "json"

def extract_json (json)
	output_array = JSON.parse(File.open(json).read)
	return output_array
end

frags_original_order = extract_json('frag_ids_original_order.json')
frags_by_density = extract_json('frags_by_density.json')
id_pos_hash = extract_json('pos_hash.json') #remember the fragments are in a random order here

