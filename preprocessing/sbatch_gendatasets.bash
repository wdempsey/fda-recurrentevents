#!/bin/bash
echo "Generating dataset"
sbatch -o MCsim.stdout.txt \
       -e MCsim.error.stdout.txt \
       --job-name=gendata \
       gendatasets_loop.sbatch


