<<'COMMENT'

This script will run algorithm_test.rb 40 times with varying parameters (4 parameter groups).
ARGV Parameters: dataset, run, gen, pop_size, select_num, c_mut, s_mut, save, ran (each explained in algorithm_test.rb).
Each for loop contains runs of the algorithm with a different div value

COMMENT

source R-3.1.0
source ruby-2.0.0

for i in {1..10}
do
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4b p_run$i 10000000000 100 50 35 35 25 5"
	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4b p_run$(($i+10)) 10000000000 20 10 7 7 5 1"
	bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4b p_run$(($i+20)) 10000000000 50 25 20 20 6 4"
	bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4b p_run$(($i+30)) 10000000000 20 4 9 9 1 1"
done