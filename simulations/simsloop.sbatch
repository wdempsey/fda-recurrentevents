#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=1000m
#SBATCH --time= 00-1:00:00              # Runtime in D-HH:MM
#SBATCH --output=/home/wdem/fda-recurrentevents/simulations/%u/%x-%j.log

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Simple R setup
module load R/3.6.3-fasrc01

# Run program
Rscript basic_simulation.R
