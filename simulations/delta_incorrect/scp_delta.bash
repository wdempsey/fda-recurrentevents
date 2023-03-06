#!/bin/bash
for ITER in {1..2}; do
    #
    echo "${ITER}"
    #
    sbatch --job-name=my_analysis_${ITER} sim_delta_simsloop.sbatch
    #
    sleep 5 # pause to be kind to the scheduler
done
