# Input parameters
# The number of asteroids processed in each Python program, e.g. 10000
batch_size=10000
# The number of batches run in parallel in each large job, e.g. 40
num_batch=25
# The smallest asteroid number to process
# min_ast_num=0
min_ast_num=1000000
# The largest asteroid number to process
# max_ast_num=542000
max_ast_num=1255514

# Compute number of jobs required to process all data
ast_per_job=$((num_batch*batch_size))
job_num_min=$((min_ast_num/ast_per_job + 1))
job_num_max=$((max_ast_num/ast_per_job + 1))
num_jobs=$((job_num_max - job_num_min + 1))
echo "Running $num_jobs jobs with batch_size=$batch_size and $num_batch batches per job..."

for ((job_num=job_num_min; job_num<=job_num_max; job_num++))
do
	bash scripts/ztf_nearest_ast_batch.sh $batch_size $num_batch $max_ast_num $job_num
	# echo "job_num=$job_num"
done
