#encoding: utf-8
class WriteIt

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
		IO.foreach(file).collect {|l| l.chomp }
	end

	# Input: file location (file must have integers on each line)
	# Output: array with each line of the file as an entry, converted to integer
	def self.file_to_ints_array(file)
	  begin 
      WriteIt.file_to_array(file).collect {|l| Integer(l) }
    resuce ArgumentError
      $stderr.puts "Not all lines in file can be converted to ints"
      exit
    end
	end

	# Input: file location (file must have floats on each line)
	# Output: array with each line of the file as an entry, converted to float
	def self.file_to_floats_array(file)
    begin 
      WriteIt.file_to_array(file).collect {|l| Float(l)}
    rescue ArgumentError
      $stderr.puts "Not all lines in file can be converted to ints"
    end
	end
end

