#!/bin/bash

for $BS in {1..2}
do
  for $MI in {1..2}
  do
    for $DELTA in 5 15 30
    do
      sbatch bootstrap_estimation_loop.sbatch
    done
  done
done
