#!/bin/bash
for ITER in {1..2}; do
    #
    echo "${ITER}"
    #
    sbatch --job-name=my_analysis_${ITER} simsloop.sbatch
    #
    sleep 10 # pause to be kind to the scheduler
done
