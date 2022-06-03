#!/bin/bash

# This script is used to submit the jobs of MEPS.R to cluster.

#$ -M zwu56@emory.edu
#$ -o MEPS.log
#$ -N MEPS
#$ -e MEPS.err
#$ -q short.q

R --vanilla  < /home/zwu56/MEPS/MEPS.R
