#!/bin/bash -xuv
##########
# Creates a CONDOR submission script to processes one configuration of one snapshot of MXXL
##########

simdir=`readlink -f $1`
outdir=`readlink -f $2`

codedir=`pwd`

condorjob=$outdir/condor.submit

#nfiles=`ls $simdir/halo_cid*.convergence_map | wc -l`
nfiles=`cat $simdir/numfiles`


if [ -e $condorjob ]; then
    rm $condorjob
fi

echo "executable = $codedir/nfwfitwrapper.sh" >> $condorjob
echo "universe = vanilla" >> $condorjob
echo "Error = $outdir/halo."'$(Process).stderr' >> $condorjob
echo "Output = $outdir/halo."'$(Process).stdout' >> $condorjob
echo "Log = $outdir/halo."'$(Process).batch.log' >> $condorjob
echo "Arguments = $simdir/halo_cid"'$(PROCESS)'" $outdir/config.sh $outdir/halo_cid"'$(PROCESS).out' >> $condorjob
echo "queue $nfiles" >> $condorjob



