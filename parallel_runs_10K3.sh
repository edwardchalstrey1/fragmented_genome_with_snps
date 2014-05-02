### ARGV Parameters ###
### dataset, run, gen, pop_size, select_num, mut_num, save, ran, div

ruby algorithm_test.rb 10K_dataset3 p_run1 10000000000 100 50 10 10 5 1000.0 & # TODO change to 100 as 1000 doesn't work, but check with small first
ruby algorithm_test.rb 10K_dataset3 p_run2 10000000000 20 10 2 2 1 100.0 &
ruby algorithm_test.rb 10K_dataset3 p_run2 10000000000 20 10 2 2 1 1000.0 