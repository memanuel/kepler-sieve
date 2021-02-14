# Start timer
t0=$SECONDS

sleep 1

# End timer
t1=$SECONDS
et_sec=$(($t1-$t0))

elapsed_time() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 mf=$(bc <<< "${1}/60")
 printf "Elapsed Time: %02d:%02d:%02d (%05.2f minutes)\n" $h $m $s
}
echo $(elapsed_time et_sec)

# et_min=$(($et_sec/60))
# et_sec=$(($et_sec-$et_min*60))
# echo "Elapsed time $et_min:$et_sec"
# date -d@36 -u +%H:%M:%S