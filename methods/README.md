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
3. 
4. The `plot_generator.R` [file](/plot_generator.R) generates plots
using the calculated numerators and denominators.
5. Save directory for plots is `/visualization/figs/`
