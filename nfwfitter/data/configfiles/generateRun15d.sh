#!/bin/bash -u


sampling=0_0
binning=hstnoisebins
r=5


for snap in 124 141; do

    cat ../shearprofiles/reduced_coresizeindex.list | { while read cluster coresizeindex redshift; do 

	    for mc in c4 duffy; do

		for xray in NONE WTG; do
	    
		    config=hstnoisebins-${mc}-r${r}-xray${xray}-${cluster}

		    dir=../../bk11_lensing/snap$snap/intlength400/$config
		    realdir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		    mkdir -p $realdir
		    ln -s $realdir $dir

		    cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
		    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
		    echo "xraycentering=${xray}" >> $dir/config.sh
	    
		    echo $config >> ../run15d.$snap

		done
		



		config=hstnoisebins-${mc}-r${r}-core${coresizeindex}-${cluster}

		dir=../../bk11_lensing/snap$snap/intlength400/$config
		realdir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		mkdir -p $realdir
		ln -s $realdir $dir

		cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh core_${coresizeindex}.sh ${binning}.sh > $dir/config.sh
		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
	    	echo "targetz=${redshift}" >> $dir/config.sh
		echo $config >> ../run15d.$snap


	    done

	done
    }
    
done

	



	
