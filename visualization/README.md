# Data Visualization #

This directory contains all the files you need to visualize the data. These methods are: 

* Parallelized script to build numerator and denominator for kernel mean estimator per participant 
* Bash script to submit jobs to server per participant
* Aggregating script to build per participant and average across participant mean curves.

## Code description ##
Steps to run code are as follows:
1. [Bash script](/sbatch_plots_loop.bash) submits a sequence of jobs
   (one for each participant).
  * Log in to the Harvard RC system
  * Find the `fda-recurrentevents` directory
  * Enter in command line: ```bash sbatch_plots_loop.bash```
2. Each job uses
   [`weight_calculator_per_participant.R`](/weight_calculator_per_participant.R)
   to calculate numerator and denominator weights at event and 1000
   sampled non-event times.
3. Save directory is `/visualization/data/`
4. The `plot_generator.R` [file](/plot_generator.R) generates plots
using the calculated numerators and denominators.
5. Save directory for plots is `/visualization/figs/`
