# Data Visualization #

This directory contains all the files you need to implement the
subsampling methodology. These methods are: 

## Code description ##
Steps to run code are as follows:
1. [Bash script](/sbatch_gendatasets.bash) submits a single job
  * Log in to the Harvard RC system
  * Find the `fda-recurrentevents` directory
  * Enter in command line: ```bash sbatch_gendatasets.bash```
2. The job runs [`generate_dataset.R`](/generate_dataset.R) to
   generate the dataset at event and a random subsampled non-event
   times. These files are saved in `/data-for-fda/data/` folder as 
3. [Mean and covariance estimation file](/mean_cov_estimation.R) can
   be run to generate `mean_estimates` and `pooled_estimates`
