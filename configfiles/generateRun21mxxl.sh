#!/bin/bash -u

snap=54
binning=hstnoisebins


cat hst_siminput.list | awk '($1 !~ "#"){print}' | { 

    while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
	
	for r in 5 16; do
	    
	    for center in xraymag szanalytic; do
	    
		for mc in c4 duffy diemer15; do

		    config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    
		
		    echo $config >> ../run21mxxl$snap
		
		    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		    
		    
		    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh core_${coresizeindex}.sh > $dir/config.sh
		    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
		    echo "targetz=$zcluster" >> $dir/config.sh
		    
		    if [ "$center" = "xraymag" ]; then
			echo "xraycentering=magneticum" >> $dir/config.sh
		    elif [ "$center" = "szanalytic" ]; then
			echo "sztheoreticalcentering=analytic" >> $dir/config.sh
		    fi
		

		done
   		
	    done
	    
	done
	
    done
    
}





	
