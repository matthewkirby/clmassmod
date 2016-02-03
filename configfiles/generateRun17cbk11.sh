#!/bin/bash -u

snap=124
binning=hstnoisebins

cat ../shearprofiles/coresizeindex.list | { 

    while read cluster coresizeindex zcluster; do 
	
	for r in 5 16; do
	    
	    for mc in c4 duffy diemer15; do

		center=xraylensingvoigt
		
		config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    
		
		echo $config >> ../run17cbk11$snap
		
		
		dir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		
		if [ ! -e $dir ]; then
		    mkdir $dir
		fi
		    
		    
		    
		cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
		echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
		echo "targetz=$zcluster" >> $dir/config.sh
		
		echo "xraycentering=voigtlensingpeak" >> $dir/config.sh
		    
   		
	    done
	    
	done
	
    done
    
}





	
