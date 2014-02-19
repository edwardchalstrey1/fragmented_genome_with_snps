#!/usr/bin/ruby
# encoding: utf-8
#
#  remove_from_karyotype.rb
#
#  Created by Dan MacLean (TSL) on 2013-06-19.
#  Copyright (c). All rights reserved.
#
require 'pp'
require 'json'
require 'csv'
require 'bio'
require 'barmcakes'

keepers = {}
File.open(ARGV[0]).each_line do |x|
  #$stderr.puts x.chomp
  keepers[x.chomp] = 1
end

File.open(ARGV[1]).each_line do |line|
 #$stderr.puts line.split(/\s/)[2]
 #$stderr.puts "yup" if keepers.has_key?(line.split(/\s/)[2])
  puts line if keepers.has_key?(line.split(/\s/)[2])
end
