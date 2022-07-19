#!/bin/bash

#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --mail-user=wdem@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000m 
#SBATCH --time=1-00:00              # Runtime in D-HH:MM
#SBATCH --account=wdem
#SBATCH --partition=standard

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Simple R setup
module load R/3.3.3-fasrc01

# Run program
Rscript bootstrap_mi.R $BOOTSTRAP
