#!/bin/bash

#SBATCH --mail-user=wdem@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=01:00:00
#SBATCH --array=1-500
#SBATCH --account=wdem1
#SBATCH --partition=standard
#SBATCH --output=/home/wdem/fda-recurrentevents/simulations/delta_incorrect/attempt2.log

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Simple R setup
module load R/4.2.0

# Run program
Rscript --verbose sim_delta.R 36 40 44 48 52 > ./outfiles/sim_delta_v2_$SLURM_ARRAY_TASK_ID.Rout
