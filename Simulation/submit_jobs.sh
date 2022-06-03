#!/bin/bash

# This script can take two arguments, SCRIPT and ANALYSIS.

##################### Change the constants as needed ##############################
analysis=$2      # change for every analysis you run (2nd arg)
maildom='@emory.edu'   # your email domain (for receiving error messages)
myscratch="./scratch"  # location of your persistent scratch dir
resultdir="./out"  # This is a folder in permanent storage
script=$1      # your code as (R or Python) script (1st arg)
max_jobs=400   # max number of jobs to run at a time
total_jobs=32000 # needs to be divisible by max_jobs
################### typically don't have to change anything below here ##############

username=$(id -nu)

# run R script to get list size
# this apparently doesn't work
# listsize=$(${script} listsize)

nloops=$((${total_jobs}/${max_jobs}-1))

# if scratch directory doesn't exist, make it
[ ! -d ${myscratch} ] && mkdir ${myscratch}
[ ! -d ${myscratch}/out ] && mkdir ${myscratch}/out
[ ! -d ${myscratch}/err ] && mkdir ${myscratch}/err

# submit first batch of jobs
for i in $(seq 1 $max_jobs); do
	qsub -N ${analysis}$i \
	-v job_control='run' \
	-v iter=$i \
	-e ${myscratch}/err \
	-o ${myscratch}/out \
	${script}
done

# sub next batches of jobs holding 
for j in $(seq 1 $nloops); do
	for i in $(seq 1 $max_jobs); do
		hold_name=${analysis}$((($j-1)*$max_jobs+$i))
		jid=$(($j*$max_jobs+$i))
		qsub -N ${analysis}$jid \
		-v job_control='run' \
		-v iter=$jid \
		-e ${myscratch}/err \
		-o ${myscratch}/out \
		-hold_jid $hold_name
		${script}
	done
done

# sub job that does merging
set -euo pipefail
last_batch=$(
for number in $(seq $(($total_jobs - $max_jobs)) $total_jobs); do
    printf '%s\n' "${analysis}${number}"
done | paste -sd ',' -
)

# want this one to wait until every one of the last jobs finishes merging
qsub -N ${analysis}_merge \
	 -v job_control='merge' \
	 -e ${myscratch}/err \
	 -o ${myscratch}/out \
	 -hold_jid $last_batch \
	 -m ea -M "${username}${maildom}" \
	 ${script}





