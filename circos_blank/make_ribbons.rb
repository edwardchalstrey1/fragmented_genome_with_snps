#!/usr/bin/ruby
# encoding: utf-8
#
#  make_ribbons.rb
#
#  Created by Dan MacLean (TSL) on 2013-06-18.
#  Copyright (c). All rights reserved.
#
require 'pp'
require 'json'
require 'csv'
require 'bio'
require 'barmcakes'

name_mapping = {"gi|436562701|gb|KB205939.1|" => "as1",
"gi|436562697|gb|KB205940.1|" => "as2", 
"gi|436562693|gb|KB205941.1|" => "as3",
"gi|436562683|gb|KB205942.1|" => "as4", 
"gi|436562679|gb|KB205943.1|" => "as5", 
"gi|436562675|gb|KB205944.1|" => "as6", 
"gi|436562672|gb|KB205945.1|" => "as6", 
"gi|436562666|gb|KB205946.1|" => "as7", 
"gi|436562663|gb|KB205947.1|" => "as8", 
"gi|436562659|gb|KB205948.1|" => "as9",
"gi|436562655|gb|KB205949.1|" => "as10",
"gi|436562653|gb|KB205950.1|" => "as11", 
"gi|436562648|gb|KB205951.1|" => "as12"
}
lines = []

Dir.glob("#{ARGV[0]}/Map_*").each do |dir|
  
File.open("#{dir}/AlignDetails.tab").each_with_index do |line,index|
  tmp = line.chomp.split(/\t/)
  id = "homology_seg" + ( lines.length / 2).to_s
  lines << "#{id} #{name_mapping[tmp[0]]} #{tmp[1]} #{tmp[2]}"
  lines << "#{id} #{tmp[3].gsub(/Cf746836_TGAC_s1v1_scaffold_/,"cf")} #{tmp[4]} #{tmp[5]}"
  
end

end
puts lines