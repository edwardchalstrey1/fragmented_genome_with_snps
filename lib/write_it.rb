#encoding: utf-8
class WriteIt
	# Input 0: Filename by which to save an array to a .txt file, one value per line
	# Input 1: Array to save
	def self.write_txt (filename, array)
		File.open("#{filename}.txt", "w+") do |f|
			array.each { |i| f.puts(i) }
		end
	end

	# Input 0: Filename by which to save an array with filetype extension, one value per line
	# Input 1: Array to save
	def self.write_data (filename, array)
		File.open("#{filename}", "w+") do |f|
			array.each { |i| f.puts(i) }
		end
	end

	# Input: An integer you wish to know whether it's a prime
	# Output: true/false
	def self.prime? (n)
		for d in 2..(n - 1)
			if (n % d) == 0
	    		return false
	    	end
		end
		true
	end

end

