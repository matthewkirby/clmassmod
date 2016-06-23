#!/bin/bash -u
#############

jobdir=$1
torun=$2
snap=$3
delta=$4
nbins=$5


#chainbase=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400
outputdir=/vol/euclid1/euclid1_2/dapple/rundlns/bk11snap$snap

jobfile=$jobdir/rundln$snap.$delta.condor.submit

if [ -e $jobfile ]; then
    rm $jobfile
fi

startbin=0
if [ "$nbins" -eq "-1" ]; then
    startbin=-1
fi


for config in `cat $torun`; do

    chaindir=`grep "bk11snap$snap $config" haloprocessed | awk '{print $3}'`

    if [ -z "$chaindir" ]; then
	echo "Cannot find $config"
	continue
    fi
    
#    chaindir=$chaindir/$config
    workdir=$outputdir/$config
    mkdir -p $workdir

    for massbin in `seq $startbin $nbins`; do
    
    
	outfile=$workdir/rundln$snap.$delta.$massbin
    
	sed -e "s|JOBDIR|$jobdir|g" -e "s|CHAINDIR|$chaindir|g" -e "s|OUTFILE|$outfile|g" -e "s|CONFIG|$config|g" -e "s|SNAP|$snap|g" -e "s|MASSBIN|$massbin|g" -e "s|DELTA|$delta|g" rundln.condor.bk11.template >> $jobfile

    done
    
done




