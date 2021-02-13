#!/bin/bash
# *****************************************************************************
# Integrate one multiple batches of asteroids
# Call from kepler-sieve/src
# $ scripts/asteroid_traj_batch min_ast_num max_ast_num
# $ scripts/asteroid_traj_batch 0 550000
# When called with job_num 1, this will integrate asteroids 0 to 1000
# With job_num 2, it will integrate 1000 to 2000, etc.

# *****************************************************************************
# Input parameters
# The smallest asteroid ID to process, e.g. 0 or 1000000
min_ast_num=$1
# The largest asteroid ID to process, e.g. 550000 or 1500000
max_ast_num=$2
# The number of asteroids processed in each Python program, e.g. 1000
batch_size=1000
# The number of batches run in parallel in each large job, e.g. 40
num_batch=40

# *****************************************************************************
# Compute number of jobs required to process all data
ast_per_job=$((num_batch*batch_size))
job_num_min=$((min_ast_num/ast_per_job + 1))
job_num_max=$((max_ast_num/ast_per_job + 1))
echo "Running $num_jobs jobs with batch_size=$batch_size and $num_batch batches per job..."

# Delegate to asteroid_traj_batch.sh
for ((job_num=job_num_min; job_num<=job_num_max; job_num++))
do
	bash scripts/asteroid_traj_batch.sh $batch_size $num_batch $max_ast_num $job_num
	# echo "job_num=$job_num"
done
