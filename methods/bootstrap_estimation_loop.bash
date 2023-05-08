#!/bin/bash

for BS in {1..2}
do
  for MI in 1
  do
    for DELTA in 15
    do
      echo 'This is BS number $BS with MI $MI and Delta equals $DELTA'
      export BS=$BS,MI=$MI,DELTA=$DELTA
      #sbatch --export=BS=$BS,MI=$MI,DELTA=$DELTA bootstrap_estimation_loop.sbatch
    done
  done
done

