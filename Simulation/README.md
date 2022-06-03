__README for simulation code__

# R package 

The `SuperLearner` package can be installed from CRAN.


# Overview 

The code for this simulation were run on a 
Linux system with a Slurm Workload Manager. The scripts can be 
executed with different commands passed in from the system 
environment. We have also included a shell script (`submit_jobs.sh`) that shows the 
workflow of executing the R scripts. However, note that some of 
the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. To read more about this work flow see 
[this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 

This folder includes the script for simulations `simulation.R`, functions
for two-stage super learner `twostageSL.R`, functions for candidate estimators `Estimators.R`
and shell script for executing simulation scripts `submit_jobs.sh`. The basic idea is that 
in the bash script we first assign the maximum number of jobs to run 
at a time (`"max_jobs"`) and the total number of jobs (`"total_jobs"`). Then we submit a 
job that computes how many loops of running will be needed.
Next `simulation.R` submitted the first batch of jobs with argument `"run"`, which executes the
number of jobs assigned in `"max_jobs"` parallelly. The `simulation.R` script source the 
two-stage super learner function stored in `twostageSL` and candidate estimators stored 
in `Estimators`. Then the `simulation.R` submitted the next batch of jobs holding, the next
batch of jobs will be executed automatically when the previous batch of jobs are complete.
Finally `simulation.R` is submitted with argument `"merge"`, 
which merges the results of the individual jobs, formats the output and
produces the tables included in the manuscript and supplement. 

# Questions

There is a substantial amount of code associated with this project and
there was a significant amount of computational burden in executing the
production size jobs. The code submitted along with the manuscript needs 
to be modified in some places in order to ease its use on different systems. 
I have not debugged the code across all different systems and I 
cannot guarantee that the code will run error-free on new systems. If you come 
across any issues, [please reach out to me by email](ziyue.wu@emory.edu) 
and I am happy to help sort them out. 
