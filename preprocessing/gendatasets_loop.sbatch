#!/bin/bash

#SBATCH --mail-user=wdem@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=01:00:00
#SBATCH --account=wdem1
#SBATCH --partition=standard
#SBATCH --output=/home/wdem/fda-recurrentevents/preprocessing/gendatasset.log

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Source the modules
source new-modules.sh

# Load required software modules

# Simple R setup
module load R/4.2.0

# Run program
Rscript generate_dataset.R 5 15 30 > ./gendatasets.Rout

