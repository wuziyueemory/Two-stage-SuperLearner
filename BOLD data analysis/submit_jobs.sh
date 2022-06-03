#!/bin/bash

# This script is used to submit the jobs of BOLD.R to cluster.

#$ -M zwu56@emory.edu
#$ -o BOLD.log
#$ -N BOLD
#$ -e BOLD.err
#$ -q short.q

R --vanilla  < /home/zwu56/BOLD/BOLD.R
