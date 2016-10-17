#!/bin/bash -u


sampling=0_0
binning=hstnoisebins
r=5
mc=c4

for snap in 41 54; do

    cat ../shearprofiles/coresizeindex.list | { while read cluster coresizeindex; do 

	    for xray in WTG CCCP SPTHST NONE; do
	    
		config=hstnoisebins-${mc}-r${r}-xray${xray}-${cluster}
		dir=../../mxxl_lensing/mxxlsnap$snap/$config
		mkdir $dir
		cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh

		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
		echo "xraycentering=${xray}" >> $dir/config.sh
	    
		echo $config >> ../run15a.$snap
	    


	    done

	done
    }

done


	
