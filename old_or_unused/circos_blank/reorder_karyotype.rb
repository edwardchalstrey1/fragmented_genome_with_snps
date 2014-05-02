#!/usr/bin/ruby
# encoding: utf-8
#
#  reorder_karyotype.rb
#
#  Created by Dan MacLean (TSL) on 2013-06-19.
#  Copyright (c). All rights reserved.
#
require 'pp'
require 'json'
require 'csv'
require 'bio'
require 'barmcakes'
order  = Hash.new {|h,k| h[k] = [] }

lines = {}
File.open(ARGV[0]).each do |line|
  l = line.split(/\s/)
  lines[l[2]] = line
end
#pp lines
done = []

File.open(ARGV[1]).lines.each_slice(2) do |a,b|
  a1 = a.split(/\s/)
  a2 = b.split(/\s/)
  #pp a2
  #pp a2
  print lines[a1[1]] #unless done.include?(a1[1])
  print lines[a2[1]] #unless done.include?(a2[1])
  done << a1[1]
  done << a2[1]
  
end



