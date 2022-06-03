__README for MEPS data analysis code__

# R package 

The `SuperLearner` package can be installed from CRAN.


# Overview 
In this folder we includes codes for Medical Expenditure Panel Survey (MEPS) data analysis.
MEPS is a national survey on the financing and use of medical care of families and individuals, 
their medical providers, and employers across the United States. 
We have also included a shell script (`submit_jobs.sh`) that shows the 
workflow of executing the R scripts `MEPS.R`. However, note that some of 
the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. To read more about this work flow see 
[this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 

This folder includes the training dataset `train.csv` which contains data from MEPS 2016,
testing dataset `test.csv` which contains data from MEPS 2017, R script for MEPS 
data analysis `MEPS.R`, functions for two-stage super learner `twostageSL.R`, functions 
for candidate estimators `Estimators.R` and shell script for executing MEPS data analysis 
R scripts `submit_jobs.sh`.
