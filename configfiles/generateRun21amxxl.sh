#!/bin/bash -u

snap=41
binning=hstnoisebins
center=xraymag
mc=diemer15

cat hst_siminput.list | awk '($1 !~ "#"){print}' | { 

    while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
	
	for r in 5 16; do
	    


	    config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    
		
	    echo $config >> ../run21amxxl$snap
	    
	    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
	    if [ ! -e $dir ]; then
		mkdir $dir
	    fi
		    
		    
		    
	    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh core_${coresizeindex}.sh > $dir/config.sh
	    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
	    echo "targetz=$zcluster" >> $dir/config.sh
	    
	    echo "xraycentering=magneticum" >> $dir/config.sh			
		
   	    	    
	done
	
    done
    
}





	
