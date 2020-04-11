#!/bin/bash

# Generic shell script for submitting many thousand short running jobs to a gridEngine
# cluster. The script goes through 2 phases of a typical HPC pipeline: 1. Parallel 
# execution, 2. merge of data. It assumes that you want to launch another script 
# (e.g. using R) which is set in the SCRIPT variable below. The number of parallel
# jobs run at each time is set through the max_jobs variable while the total jobs
# is set through the total_jobs variable. The first batch of jobs are executed 
# parallely while holding the second batch of jobs, when the first batch of jobs
# are done, the second batch of jobs while be automatically executed. Note the 
# total_jobs need to be divisible by max_jobs to ensure sequentially execuation. 

# This script can take two arguments, SCRIPT and ANALYSIS.

##################### Change the constants as needed ##############################
analysis=$2      # change for every analysis you run (2nd arg)
maildom='@emory.edu'   # your email domain (for receiving error messages)
myscratch="./scratch"  # location of your persistent scratch dir
resultdir="./out"  # This is a folder in permanent storage
script=$1      # your code as (R or Python) script (1st arg)
max_jobs=400   # max number of jobs to run at a time
total_jobs=3200 # needs to be divisible by max_jobs
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





