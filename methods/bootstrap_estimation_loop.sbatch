#!/bin/bash

#SBATCH --mail-user=wdem@umich.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000m
#SBATCH --time=10:00:00
#SBATCH --account=wdem1
#SBATCH --partition=standard
#SBATCH --output=/home/wdem/fda-recurrentevents/methods/attempt_estimation.log

cd $SLURM_SUBMIT_DIR

echo $SLURM_SUBMIT_DIR

echo 'Bootstrap number $BS'

echo 'MI number $MI'

echo 'Delta equals $DELTA'


# Load required software modules

# Simple R setup
module load R/4.2.0
module load gsl

# Run program
Rscript bootstrap_mi.R $BS $MI $DELTA > ./outfiles/sim_estimation_bs_$BS_mi_$MI_Delta_$DELTA.Rout
