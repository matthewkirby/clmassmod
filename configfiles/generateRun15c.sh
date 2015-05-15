#!/bin/bash -u


sampling=0_0
binning=hstnoisebins
r=5

for snap in 41 54; do

    cat ../shearprofiles/coresizeindex.list | { while read cluster coresizeindex redshift; do 

	    for mc in duffy c4; do

	    
		config=hstnoisebins-${mc}-r${r}-core${coresizeindex}-${cluster}
		dir=../../mxxl_lensing/mxxlsnap$snap/$config
		if [ -e $dir ]; then
		    mv $dir ../../mxxl_lensing/mxxlsnap$snap/old_${config}
		fi
		realdir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		mkdir $realdir
		ln -s $realdir $dir
		cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh core_${coresizeindex}.sh ${binning}.sh > $realdir/config.sh

		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $realdir/config.sh
	    	echo "targetz=${redshift}" >> $realdir/config.sh
		echo $config >> ../run15c.$snap


		config=hstnoisebins-${mc}-r${r}-sztcenter-${cluster}
		dir=../../mxxl_lensing/mxxlsnap$snap/$config
		if [ -e $dir ]; then
		    mv $dir ../../mxxl_lensing/mxxlsnap$snap/old_${config}
		fi
		realdir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		mkdir $realdir
		ln -s $realdir $dir
		cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh > $realdir/config.sh

		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $realdir/config.sh
		echo "sztheoreticalcentering=True" >> $realdir/config.sh
		echo "targetz=${redshift}" >> $realdir/config.sh
	    
		echo $config >> ../run15c.$snap

	    


	    done


	done
    }

done


	
