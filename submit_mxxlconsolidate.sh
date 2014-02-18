#!/bin/bash
#############

jobdir=$1
torun=$2

if [ -e $resubdir/consolidate.condor.submit ]; then
    rm $resubdir/consolidate.condor.submit
fi

for config in `cat $torun`; do

    sed -e "s|JOBDIR|$jobdir|g" -e "s/CONFIG/$config/g" condor.consolidate.template >> $jobdir/consolidate.condor.submit

 
done
