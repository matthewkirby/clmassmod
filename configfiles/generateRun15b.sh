#!/bin/bash -u


sampling=0_0
binning=hstnoisebins
r=5
mc=duffy

for snap in 41 54; do

    cat ../shearprofiles/reduced_coresizeindex.list | { while read cluster coresizeindex; do 

	    for xray in WTG CCCP SPTHST NONE; do
	    
		config=hstnoisebins-${mc}-r${r}-xray${xray}-${cluster}
		dir=../../mxxl_lensing/mxxlsnap$snap/$config
		if [ -e $dir ]; then
		    mv $dir ../../mxxl_lensing/mxxlsnap$snap/old_${config}
		fi
		realdir=/vol/euclid1/euclid1_3/dapple/mxxl_lensing/mxxlsnap$snap/$config
		mkdir $realdir
		ln -s $realdir $dir
		cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh > $realdir/config.sh

		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $realdir/config.sh
		echo "xraycentering=${xray}" >> $realdir/config.sh
	    
		echo $config >> ../run15b.$snap
	    


	    done

	done
    }

done


	
