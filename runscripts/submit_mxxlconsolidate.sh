#!/bin/bash
#############

jobdir=$1
torun=$2

if [ -e $jobdir/consolidate.condor.submit ]; then
    rm $jobdir/consolidate*.condor.submit
fi

for snap in 41 54; do

    for config in `cat $torun`; do

	sed -e "s|JOBDIR|$jobdir|g" -e "s/CONFIG/$config/g" -e "s/SNAP/$snap/g" condor.consolidate.template >> $jobdir/consolidate$snap.condor.submit

 
    done

done
