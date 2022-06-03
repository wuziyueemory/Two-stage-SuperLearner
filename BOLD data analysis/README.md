__README for BOLD data analysis code__

# R package 

The `SuperLearner` package can be installed from CRAN.


# Overview 
In this folder we includes codes for Back Pain Outcomes using Longitudinal Data (BOLD) data analysis.
BOLD is a large, community-based registry of patients aged 65 years and older who presented with 
primary care visits for a new episode of back pain from March 2011 to March 2013.
We have also included a shell script (`submit_jobs.sh`) that shows the 
workflow of executing the R scripts `BOLD.R`. However, note that some of 
the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. To read more about this work flow see 
[this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 

This folder includes the training dataset `train.csv` which contains data from BOLD (served
as both training and validation set through 10-fold cross-validation), 
R script for BOLD data analysis `BOLD.R`, functions for two-stage super learner `twostageSL.R`, functions 
for candidate estimators `Estimators.R` and shell script for executing BOLD data analysis 
R scripts `submit_jobs.sh`.
