#!/bin/bash
for PARTICIPANTID in {1001..1091}; do
    #
    echo "${PARTICIPANTID}"
    export PARTICIPANTID
    #
    sbatch -o MCsim_${PARTICIPANTID}.stdout.txt \
	   -e MCsim_${PARTICIPANTID}.stdout.txt \
           --job-name=analysisof_${PARTICIPANTID} \
	   plots_loop.sbatch
    #
    sleep 1 # pause to be kind to the scheduler
done



