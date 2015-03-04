#!/bin/bash -u

snap=41
sampling=0_0
binning=hstnoisebins
r=5

cat ../shearprofiles/coresizeindex.list | { while read cluster coresizeindex; do 

	for mc in c4 duffy; do


	    
	    config=hstnoisebins-${mc}-r${r}-corexray-${cluster}
	    dir=../../mxxl_lensing/mxxlsnap$snap/$config
	    mkdir $dir
	    cat scanpdf.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh core_xray.sh ${binning}.sh > $dir/config.sh
	    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
	    
	    echo $config >> ../run11a
	    


	done
   	
    done
    
}


	
