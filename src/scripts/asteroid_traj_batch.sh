#!/bin/bash
# *****************************************************************************
# Integrate one batch of asteroid trajectories on multiple threads
# Call from kepler-sieve/src
# $ scripts/asteroid_traj_batch mode batch_size num_batch max_ast_num job_num
# $ scripts/asteroid_traj_batch CSV 1000 40 550000 1
# When called with job_num 1, this will integrate asteroids 0 to 1000
# With job_num 2, it will integrate 1000 to 2000, etc.

# *****************************************************************************
# Imports
source scripts/bash_utils.sh

# *****************************************************************************
# Extract command line arguments
# The mode; one of CSV, DB, or INS
mode=$1
# This is the number of asteroids processed in each Python program, e.g. 1000
batch_size=$2
# The number of batches run in parallel in each large job, e.g. 40
num_batch=$3
# The largest asteroid number to process, e.g. 550000
max_ast_num=$4
# The number of this job, starting from 1
job_num=$5

# *****************************************************************************
# Start timer
t0=$SECONDS

# Additional parameters
sleep_time=0.1

# Set up index ranges
# j is the multiplier of the batch size; ranges from j0 to j1
j0=$((num_batch*(job_num-1)))
j1=$((j0+num_batch-1))
# n is the first asteroid number to process in each call; ranges from n0 to n1
n0=$((j0*batch_size))
n1=$(((j1+1)*batch_size))
echo "Bash asteroid_traj_batch.sh is processing asteroids from n0=$n0 to n1=$n1 with batch_size=$batch_size in $mode mode..."

# Run all the jobs jobs in parallel
for (( i=0; i<num_batch; i++))
do
	# j is the multiplier of the batch size
	j=$((j0+i))
	# n is is the asteroid number for the current python job
	n=$((j*batch_size))
	# if n is past the largest known asteroid, terminate the loop early
	if [ $n -gt $max_ast_num ]
	then
		break
	fi

	# run the last job with a progress bar, the rest silently
	if [ $i -lt $((num_batch-1)) ]
	then		
		python asteroid_integrate.py $n $batch_size $mode --quiet &	
	else
		python asteroid_integrate.py $n $batch_size $mode &
	fi

	# Slight pause so the batches will be executed in the specified order
	sleep $sleep_time
	# Save the process ID
	pids[i]=$!
	pid=$((pids[i]))
	# echo "i=$i, n=$n, pid=$pid"
done

# Wait for all outstanding jobs to be completed
for (( i=0; i<num_batch; i++))
do
	j=$((j0+i))
	n=$((j*batch_size))
	pid=$((pids[i]))
	wait $pid
	echo "Process for i=$i, n=$n, pid=$pid is complete."
done

# End timer
t1=$SECONDS
et=$(($t1-$t0))

# *****************************************************************************
# Report results
echo -e "\n********************************************************************************"
echo "$(date +"%Y-%m-%d %H:%M:%S") Done! Processed asteroid trajectories from $n0 to $n1."
elapsed_time et
echo -e "********************************************************************************\n"
