#!/bin/bash
#############

jobdir=$1
torun=$2
nbins=$3
delta=$4

for snap in 41 54; do

    jobfile=$jobdir/rundln$snap.$delta.condor.submit

    if [ -e $jobfile ]; then
	rm $jobfile
    fi


    for config in `cat $torun`; do

	for massbin in `seq 0 $nbins`; do

	    sed -e "s|JOBDIR|$jobdir|g" -e "s/CONFIG/$config/g" -e "s/SNAP/$snap/g" -e "s/MASSBIN/$massbin/g" -e "s/DELTA/$delta/g" rundln.condor.template >> $jobfile

	done
 
    done

done
