#!/bin/bash
for BOOTSTRAP in {1..2}; do
    #
    echo "${BOOTSTRAP}"
    export BOOTSTRAP
    #
    sbatch -o bootstrapsim_${BOOTSTRAP}.stdout.txt \
	   -e bootstrapsim_${BOOTSTRAP}.stdout.txt \
           --job-name=bootstrap_analysis_${BOOTSTRAP} \
	   bootstrap_loop.sbatch
    #
    sleep 1 # pause to be kind to the scheduler
done
