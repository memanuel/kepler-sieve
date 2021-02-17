# Bash script utitities
# 2021-02-17

# ******************************************************************************
# Print the elapsed time in a user friendly format
elapsed_time() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 mf=$(bc <<< "scale=2; ${1}/60")
 printf "Elapsed Time: %02d:%02d:%02d (%5.2f minutes)\n" $h $m $s $mf
}
