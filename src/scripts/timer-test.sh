# Imports
source ./bash_utils.sh

# Start timer
t0=$SECONDS

sleep 1

# End timer
t1=$SECONDS
et_sec=$(($t1-$t0))
#et=900

# print status
elapsed_time $et_sec
