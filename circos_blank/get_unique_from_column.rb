#!/usr/bin/ruby
# encoding: utf-8
#
#  get_unique_from_column.rb
#
#  Created by Dan MacLean (TSL) on 2013-06-19.
#  Copyright (c). All rights reserved.
#
require 'pp'
require 'json'
require 'csv'
require 'bio'
require 'barmcakes'


puts File.open(ARGV[0]).list_of_things_in_column(1,"\s")