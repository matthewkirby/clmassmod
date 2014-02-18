#!/bin/bash
###############


simtype=$1
fails=$2
jobdir=$3


##############

if [ $simtype == "bcc" ]; then

    simdir='/u/ki/dapple/nfs/bcc_clusters/recentered'
    simtemplate='cluster_|HID|'
    simext="hdf5"

elif [ $simtype == "bk11snap124" ]; then

    simdir='/u/ki/dapple/nfs/beckersims/snap124/intlength400'
    simtemplate='haloid|HID|_zLens0.49925070_intlength400'
    simext="fit"

elif [ $simtype == "bk11snap141" ]; then

    simdir='/u/ki/dapple/nfs/beckersims/snap141/intlength400'
    simtemplate='haloid|HID|_zLens0.24533000_intlength400'
    simext="fit"
fi
    

cat $fails | { 

    while read fail config haloid; do

	simbase=`echo $simtemplate | sed -e"s/|HID|/$haloid/"`
	simfile=$simdir/$simbase.$simext

	outdir=$simdir/$config
    

	jobname=$simtype.$config.$haloid.redo

	configfile=$outdir/config.sh
    
	outfile=$outdir/$simbase.out
	logfile=$jobdir/$simbase.log
    
	jobfile=$jobdir/p300.$jobname
    
	echo "#!/bin/bash" > $jobfile
    
	echo "bsub -q short -oo $logfile ./nfwfit.py $simfile $configfile $outfile"  >> $jobfile
    
	chmod a+x $jobfile


    done

}