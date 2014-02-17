#!/bin/bash
#############

resubdir=$1
failfile=$2

if [ ! -e $resubdir ]; then
    mkdir $resubdir
fi

if [ -e $resubdir/condor.submit ]; then
    rm $resubdir/condor.submit
fi

cat $failfile | { 

    while read fail config haloid; do


	sed -e "s|RESUBDIR|$resubdir|g" -e "s/HALOID/$haloid/g" -e "s/CONFIG/$config/g" condor.failresub.template >> $resubdir/condor.submit

 
    done
}