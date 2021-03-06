#!/bin/bash -u

snap=41
sampling=0_0
binning=hstnoisebins
r=5
mc=c4

for profile in simple highnoise lownoise ; do


    config=hstnoisebins-${mc}-r${r}-${profile}
    dir=../../mxxl_lensing/mxxlsnap$snap/$config
    mkdir $dir
    cat scanpdf.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh core_none.sh ${binning}.sh > $dir/config.sh
    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${profile}.profile" >> $dir/config.sh
	    
    echo $config >> ../run12b
	    

    
done




	
