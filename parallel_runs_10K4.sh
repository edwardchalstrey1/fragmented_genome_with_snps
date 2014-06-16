<<'COMMENT'

This script will run algorithm_test.rb 120 times with varying parameters.
ARGV Parameters: dataset, run, gen, pop_size, select_num, c_mut, s_mut, save, ran, div (each explained in algorithm_test.rb).
Each for loop contains runs of the algorithm with a different div value

COMMENT

source R-3.1.0
source ruby-2.0.0

# for i in {1..10}
# do
# 	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 1000.0 restart"
# 	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 1000.0 restart"
# 	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 1000.0 restart"
# 	bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 1000.0 restart"
# done

# for i in {41..50}
# do
# 	bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 10000.0 restart"
# 	bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 10000.0 restart"
# 	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 10000.0 restart"
# 	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 10000.0 restart"
# done

# for i in {81..90}
# do
# 	bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 100000.0 restart"
# 	bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 100000.0 restart"
# 	bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 100000.0 restart"
# 	bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 100000.0 restart"
# done

# for i in {121..130}
# do
# 	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$i 10000000000 100 50 35 35 25 5 1000000.0"
# 	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+10)) 10000000000 20 10 7 7 5 1 1000000.0"
# 	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+20)) 10000000000 50 25 20 20 6 4 1000000.0"
# 	bsub -q TSL-Test128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run$(($i+30)) 10000000000 20 4 9 9 1 1 1000000.0"
# done

bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run34 10000000000 20 4 9 9 1 1 1000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run39 10000000000 20 4 9 9 1 1 1000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run51 10000000000 20 10 7 7 5 1 10000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run35 10000000000 20 4 9 9 1 1 1000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run59 10000000000 20 10 7 7 5 1 10000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run31 10000000000 20 4 9 9 1 1 1000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run50 10000000000 100 50 35 35 25 5 10000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run36 10000000000 20 4 9 9 1 1 1000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run54 10000000000 20 10 7 7 5 1 10000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run120 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run100 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run106 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run91 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run107 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run112 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run105 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run101 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run96 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run103 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run94 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run97 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run98 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run116 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run84 10000000000 100 50 35 35 25 5 100000.0 restart"
bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run83 10000000000 100 50 35 35 25 5 100000.0 restart"
bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run81 10000000000 100 50 35 35 25 5 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run119 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run118 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run99 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Test256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run52 10000000000 20 10 7 7 5 1 10000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run114 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run93 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run113 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run115 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run104 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run111 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run117 10000000000 20 4 9 9 1 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run102 10000000000 50 25 20 20 6 4 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run95 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run92 10000000000 20 10 7 7 5 1 100000.0 restart"
bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run85 10000000000 100 50 35 35 25 5 100000.0 restart"
bsub -q TSL-Prod128 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run82 10000000000 100 50 35 35 25 5 100000.0 restart"
bsub -q TSL-Prod256 -We 10000 "ruby algorithm_test.rb 10K_dataset4 p_run108 10000000000 50 25 20 20 6 4 100000.0 restart"




