#!/bin/bash
for ITER in {1..2}; do
    #
    echo "${ITER}"
    #
    sbatch -o MCsim_${ITER}.stdout.txt \
	   -e MCsim_${ITER}.stdout.txt \
           --job-name=my_analysis_${ITER} \
	   simsloop.sbatch
    #
    sleep 10 # pause to be kind to the scheduler
done
