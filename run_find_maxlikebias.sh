#!/bin/bash -u
#############

jobdir=$1
torun=$2
snap=$3
delta=$4
nbins=$5

#chainbase=/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap$snap
outputdir=/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap$snap/


startbin=0
if [ "$nbins" -eq "-1" ]; then
    startbin=-1
fi

for config in `cat $torun`; do
	


    chaindir=`grep "mxxlsnap$snap $config" haloprocessed | awk '{print $3}'`
	
    if [ -z "$chaindir" ]; then
	echo "Cannot find $config"
	continue
    fi

    #	chaindir=$chainbase/$config
    workdir=$outputdir/$config
    mkdir -p $workdir


    for massbin in `seq $startbin $nbins`; do

	outfile=$workdir/rundln$snap.$delta.$massbin

	python ./find_maxlikebias.py mxxlsnap${snap} $chaindir $outfile $delta $massbin
	    
	
    done
 
done


