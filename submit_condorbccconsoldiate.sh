#!/bin/bash
#############

jobdir=$1
torun=$2


for config in `cat $torun`; do

    sed -e "s|JOBDIR|$jobdir|g" -e "s/CONFIG/$config/g" condor.bcc.consolidate.template >> $jobdir/consolidatebcc.condor.submit

    
done
