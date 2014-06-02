<<'COMMENT'

This script will run algorithm_test.rb 120 times with varying parameters.
ARGV Parameters: dataset, run, gen, pop_size, select_num, c_mut, s_mut, save, ran, div (each explained in algorithm_test.rb).
Each for loop contains runs of the algorithm with a different div value

COMMENT

source R-3.1.0
source ruby-2.0.0

for i in {1..10}
do
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 1000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 1000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 1000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 1000.0"
done

for i in {41..50}
do
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 10000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 10000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 10000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 10000.0"
done

for i in {81..90}
do
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 100000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 100000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 100000.0"
	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 100000.0"
done