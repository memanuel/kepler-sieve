#!/bin/bash
# *****************************************************************************
# Integrate one multiple batches of asteroids
# Call from kepler-sieve/src
# $ scripts/asteroid_traj_batch min_ast_num max_ast_num mode
# $ scripts/asteroid_traj_batch 0 550000 DB
# $ scripts/asteroid_traj_batch 1000000 1415000 CSV
# When called with job_num 1, this will integrate asteroids 0 to 1000
# With job_num 2, it will integrate 1000 to 2000, etc.

# *****************************************************************************
# Imports
source scripts/bash_utils.sh

# *****************************************************************************
# Input parameters
# The smallest asteroid ID to process, e.g. 0 or 1000000
min_ast_num=$1
# The largest asteroid ID to process, e.g. 550000 or 1500000
max_ast_num=$2
# The mode; defaults to DB
mode=${3:-DB}
# Additional options, e.g. --run_all
suffix=$4
# The number of asteroids processed in each Python program, e.g. 1000
batch_size=250
# The number of batches run in parallel in each large job, e.g. 40
num_cpu=$(nproc)
num_batch=$((num_cpu/2))

# *****************************************************************************
# Start timer
t0=$SECONDS

# Compute number of jobs required to process all data
ast_per_job=$((num_batch*batch_size))
job_num_min=$((min_ast_num/ast_per_job + 1))
job_num_max=$((max_ast_num/ast_per_job + 1))
num_jobs=$((job_num_max-job_num_min+1))
echo "Running $num_jobs jobs in $mode mode with batch_size=$batch_size and $num_batch batches per job..."

# Delegate to asteroid_traj_batch.sh
for ((job_num=job_num_min; job_num<=job_num_max; job_num++))
do
	bash scripts/asteroid_traj_batch.sh $mode $batch_size $num_batch $max_ast_num $job_num $suffix
	# echo "job_num=$job_num"
done

# End timer
t1=$SECONDS
et=$(($t1-$t0))

# *****************************************************************************
# Report results
echo -e "\n********************************************************************************"
echo "$(date +"%Y-%m-%d %H:%M:%S") Done! Completed $num_jobs."
elapsed_time et
echo -e "********************************************************************************\n"
