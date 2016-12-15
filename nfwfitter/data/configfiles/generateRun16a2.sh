#!/bin/bash -u


sampling=0_0
binning=hstnoisebins

snap=124

center=xrayNONE

cat ../shearprofiles/coresizeindex.list | { 
    
    while read cluster coresizeindex redshift; do 

	core=core$coresizeindex
	
	for mc in c4 duffy; do

	    for r in 5 16; do
	    

	    
		config=hstnoisebins-${mc}-r${r}-${center}-${cluster}

		dir=../../bk11_lensing/snap$snap/intlength400/$config
		realdir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		mkdir -p $realdir
		ln -s $realdir $dir

		cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh

		echo $config >> ../run16a2.$snap


	    done

	done

    done
}



	



	
