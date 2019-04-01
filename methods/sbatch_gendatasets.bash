#!/bin/bash
#for PARTICIPANTID in {1001..1091}; do
    #
echo "Generating dataset"
#
sbatch -o MCsim.stdout.txt \
       -e MCsim.error.stdout.txt \
       --job-name=gendata \
       gendatasets_loop.sbatch
#
#sleep 1 # pause to be kind to the scheduler
#done



