#!/bin/bash

for BS in 1
do
    for MI in 1
    do
	for DELTA in 15
	do
            echo "Bootstrap number $BS, MI number $MI, current DELTA is $DELTA"
	    export BS=$BS,MI=$MI,DELTA=$DELTA
	    sbatch --export=BS=$BS,MI=$MI,DELTA=$DELTA bootstrap_estimation_loop.sbatch
	    sleep 5
	done
    done
done

