require "rubygems"
require "rinruby"

def scatta (dataset, method_name_string, filename_string)
	myr = RinRuby.new(echo = false)
	file = "~/fragmented_genome_with_snps/arabidopsis_datasets/"+dataset.to_s+"/re_files/id_n_density_txts/"+filename_string+".txt"
	myr.assign "file", file
	myr.eval 'd <- as.vector(as.matrix(read.table(file, quote="\"")))'
	myr.eval 'source("~/fragmented_genome_with_snps/scatter_vectors.R")'
	myr.assign "filename_string", filename_string
	myr.assign "method_name_string", method_name_string
	myr.assign "dataset", dataset.to_s
	myr.eval 'scatta(d, method_name_string, filename_string, dataset)'
end
Dir.mkdir(File.join(Dir.home, "fragmented_genome_with_snps/arabidopsis_datasets/"+ARGV[0].to_s+"/figures"))
scatta(ARGV[0], "Original
       Order", "d_o")
scatta(ARGV[0], "Control 1", "d_c1")
scatta(ARGV[0], "Control 2", "d_c2")
scatta(ARGV[0], "Method 1a", "d_m1a")
scatta(ARGV[0], "Method 1b", "d_m1b")
scatta(ARGV[0], "Method 2a 
       (left)", "d_m2al")
scatta(ARGV[0], "Method 2a 
       (right)", "d_m2ar")
scatta(ARGV[0], "Method 2a", "d_m2a")
scatta(ARGV[0], "Method 2b 
       (left)", "d_m2bl")
scatta(ARGV[0], "Method 2b 
       (right)", "d_m2br")
scatta(ARGV[0], "Method 2b", "d_m2b")
