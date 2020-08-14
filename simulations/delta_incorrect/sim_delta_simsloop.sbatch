#!/bin/bash

#SBATCH --mail-user=wdem@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=01:00:00
#SBATCH --account=wdem1
#SBATCH --partition=standard
#SBATCH --output=/home/wdem/fda-recurrentevents/simulations/simdelta/%u/%x-%j.log

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Simple R setup
module load R/3.6.1

# Run program
Rscript sim_delta.R 52
Rscript sim_delta.R 48
Rscript sim_delta.R
Rscript sim_delta.R 40
Rscript sim_delta.R 36
