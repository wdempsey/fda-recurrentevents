#!/bin/bash
for BOOTSTRAP in {1..200}; do
    #
    echo "${BOOTSTRAP}"
    export BOOTSTRAP
    #
    sbatch -o MCsim_${BOOTSTRAP}.stdout.txt \
	   -e MCsim_${BOOTSTRAP}.stdout.txt \
           --job-name=my_analysis_${BOOTSTRAP} \
	   simsloop.sbatch
    #
    sleep 1 # pause to be kind to the scheduler
done
