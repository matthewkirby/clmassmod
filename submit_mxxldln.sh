#!/bin/bash -u
#############

jobdir=$1
torun=$2
snap=$3
delta=$4
nbins=$5

chainbase=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap
outputdir=/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap$snap/

jobfile=$jobdir/rundln$snap.$delta.condor.submit

    if [ -e $jobfile ]; then
	rm $jobfile
    fi

    for config in `cat $torun`; do

	chaindir=$chainbase/$config
	workdir=$outputdir/$config
	mkdir -p $workdir


	for massbin in `seq 0 $nbins`; do

	    outfile=$workdir/rundln$snap.$delta.$massbin

	    sed -e "s|JOBDIR|$jobdir|g" -e "s|CHAINDIR|$chaindir|g" -e "s|OUTFILE|$outfile|g" -e "s|CONFIG|$config|g" -e "s|SNAP|$snap|g" -e "s|MASSBIN|$massbin|g" -e "s|DELTA|$delta|g" rundln.condor.template >> $jobfile

	done
 
    done


