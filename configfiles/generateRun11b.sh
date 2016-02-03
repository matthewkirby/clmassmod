#!/bin/bash -u


sampling=0_0
binning=hstnoisebins
r=5

for snap in 41 54; do

cat ../shearprofiles/coresizeindex.list | { while read cluster coresizeindex; do 

	for mc in c4 duffy; do


	    
	    config=hstnoisebins-${mc}-r${r}-corexraycccp-${cluster}
	    dir=../../mxxl_lensing/mxxlsnap$snap/$config
	    mkdir $dir
	    cat scanpdf.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh core_xraycccp.sh ${binning}.sh > $dir/config.sh
	    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
	    
	    echo $config >> ../run11b.$snap
	    


	done
   	
    done
    
}

done


	
