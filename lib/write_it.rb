#encoding: utf-8
class WriteIt

	require 'text-table'

	# Input 0: Filename by which to save an array to a .txt file, one value per line
	# Input 1: Array to save
	def self.write_txt(filename, array)
		File.open("#{filename}.txt", "w+") do |f|
			array.each { |i| f.puts(i) }
		end
	end

	# Input 0: Filename by which to save an array with filetype extension, one value per line
	# Input 1: Array to save
	def self.write_data(filename, array)
		File.open("#{filename}", "w+") do |f|
			array.each { |i| f.puts(i) }
		end
	end

	# Input: file location
	# Output: array with each line of the file as an entry
	def self.file_to_array(file)
		array = []
		IO.foreach(file) { |line| array << line.gsub(/\n/,'') }
		array
	end

	# Input: file location (file must have integers on each line)
	# Output: array with each line of the file as an entry, converted to integer
	def self.file_to_ints_array(file)
		array = []
		IO.foreach(file) { |line| array << line.gsub(/\n/,'').to_i }
		array
	end

	# Input: file location (file must have floats on each line)
	# Output: array with each line of the file as an entry, converted to float
	def self.file_to_floats_array(file)
		array = []
		IO.foreach(file) { |line| array << line.gsub(/\n/,'').to_f }
		array
	end
end

