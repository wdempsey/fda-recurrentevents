#!/bin/bash

#SBATCH --mail-user=wdem@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000m
#SBATCH --time=10:00:00
#SBATCH --array=1-200
#SBATCH --account=wdem1
#SBATCH --partition=standard
#SBATCH --output=/home/wdem/fda-recurrentevents/methods/attempt.log

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

# Load required software modules

# Simple R setup
module load R/4.2.0
module load gsl

# Run program
Rscript bootstrap_mi.R $SLURM_ARRAY_TASK_ID > ./outfiles/sim_$SLURM_ARRAY_TASK_ID.Rout
