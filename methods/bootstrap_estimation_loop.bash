#!/bin/bash

for BS in {1..200}
do
    for MI in {1..2}
    do
	for DELTA in 5 15 30
	do
            echo "Bootstrap number $BS, MI number $MI, current DELTA is $DELTA"
	    export BS=$BS,MI=$MI,DELTA=$DELTA
	    sbatch --export=BS=$BS,MI=$MI,DELTA=$DELTA bootstrap_estimation_loop.sbatch
	    sleep 1
	done
    done
done

