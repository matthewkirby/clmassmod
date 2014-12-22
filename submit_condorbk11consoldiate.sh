#!/bin/bash
#############

jobdir=$1
torun=$2
tag=$3

if [ -e $jobdir/consolidate.${tag}.condor.submit ]; then
    rm $jobdir/consolidate*${tag}.condor.submit
fi

for snap in 141 124; do

    for config in `cat $torun`; do

	sed -e "s|JOBDIR|$jobdir|g" -e "s/CONFIG/$config/g" -e "s/SNAP/$snap/g" condor.bk11.consolidate.template >> $jobdir/consolidate$snap.${tag}.condor.submit

 
    done

done
