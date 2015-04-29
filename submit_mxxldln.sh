#!/bin/bash
#############

jobdir=$1
torun=$2
nbins=$3
delta=$4

for snap in 41 54; do

    if [ -e $jobdir/rundln$snap.condor.submit ]; then
	rm $jobdir/rundln$snap.condor.submit
    fi


    for config in `cat $torun`; do

	for massbin in `seq 0 $nbins`; do

	    sed -e "s|JOBDIR|$jobdir|g" -e "s/CONFIG/$config/g" -e "s/SNAP/$snap/g" -e "s/MASSBIN/$massbin/g" -e "s/DELTA/$delta/g" rundln.condor.template >> $jobdir/rundln$snap.condor.submit

	done
 
    done

done
