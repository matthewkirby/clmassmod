#!/bin/bash -u

snap=41
binning=hstnoisebins

cat hst_siminput.list | awk '($1 !~ "#"){print}' | { 

    while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
	
	for r in 5 16; do
	    
	    for mc in c4 duffy diemer15; do

		for center in szanalytic szxvpbcg; do
		
		    config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    
		
		    echo $config >> ../run17dmxxl$snap
		
		    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		    
		    
		    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh core_${coresizeindex}.sh > $dir/config.sh
		    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
		    echo "targetz=$zcluster" >> $dir/config.sh

		    if [ "$center" == "szanalytic" ]; then
			echo "sztheoreticalcentering=analytic" >> $dir/config.sh
			echo "szbeam=1.2" >> $dir/config.sh
			echo "sz_xi=$spt_xi" >> $dir/config.sh
		    elif [ "$center" == "szxvpbcg" ]; then
			echo "sztheoreticalcentering=xvp_szbcg" >> $dir/config.sh
		    fi
			
		

		done
   		
	    done
	    
	done
	
    done
    
}





	
