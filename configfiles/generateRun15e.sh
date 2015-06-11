#!/bin/bash -u


sampling=0_0
binning=hstnoisebins
r=5

snap=124



cat ../shearprofiles/reduced_coresizeindex.list | { 
    
    while read cluster coresizeindex redshift; do 
	
	for mc in c4 duffy; do
	    
	    for center in sztcenter xraySPTHST; do
	    
		config=hstnoisebins-${mc}-r${r}-${center}-${cluster}

		dir=../../bk11_lensing/snap$snap/intlength400/$config
		realdir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		mkdir -p $realdir
		ln -s $realdir $dir

		cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh

		if [ ${center} == "xraySPTHST" ]; then
		    echo "xraycentering=SPTHST" >> $dir/config.sh
		else
		    echo "sztheoreticalcentering=True" >> $dir/config.sh
	    	    echo "targetz=${redshift}" >> $dir/config.sh
		fi


	    
		echo $config >> ../run15e.$snap

	    done

	done

    done
}



	



	
