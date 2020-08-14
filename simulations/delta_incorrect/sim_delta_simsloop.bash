#!/bin/bash
for ITER in {1..1000}; do
    #
    echo "${ITER}"
    #
    sbatch --job-name=my_analysis_${ITER} simsloop.sbatch
    #
    sleep 5 # pause to be kind to the scheduler
done
