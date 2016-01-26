#!/bin/bash -u


binning=hstnoisebins
r=5

for snap in 41 54; do

    cat hst_siminput.list | awk '($1 !~ "#"){print}' | { 

	while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
		    
	    for center in xrayNONE xraymag szanalytic szmag; do
		
		for mc in duffy c4 diemer15; do
		    
		    config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    
		
		    echo $config >> ../run22mxxl$snap
		
		    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		   
		 
		    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh core_${coresizeindex}.sh > $dir/config.sh

		    echo "targetz=$zcluster" >> $dir/config.sh
		    echo "szbeam=1.2" >> $dir/config.sh
		    echo "sz_xi=$spt_xi" >> $dir/config.sh


		    if [ "$center" = "xrayNONE" ]; then
			echo "xraycentering=None" >> $dir/config.sh
			echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh		
		    elif [ "$center" = "xraymag" ]; then
			echo "xraycentering=magneticum" >> $dir/config.sh
			echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.xraycenter.profile" >> $dir/config.sh		
		    elif [ "$center" = "szanalytic" ]; then
			echo "sztheoreticalcentering=analytic" >> $dir/config.sh
			echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh		
		    elif [ "$center" = "szmag" ]; then
			echo "sztheoreticalcentering=magneticum" >> $dir/config.sh
			echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh		
		    fi


		done
   		
	    done
	   	
	done
    
    }

done


binning=lineargaussbins12
shapenoise=0.25
r=9

for snap in 41 54; do

    cat megacam_siminput.list | awk '($1 !~ "#"){print}' | { 

	while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
		    
	    for center in xrayNONE xraymag szanalytic szmag; do
		
		for mc in duffy c4 diemer15; do

		    config=mega-${mc}-r${r}-sigma${shapenoise}-${center}-${cluster}	    
		
		    echo $config >> ../run22mxxl$snap
		
		    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		   
		 
		    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh core_${coresizeindex}.sh > $dir/config.sh

		    echo "targetz=$zcluster" >> $dir/config.sh
		    echo "nperarcmin=$ndensity" >> $dir/config.sh
		    echo "shapenoise=$shapenoise" >> $dir/config.sh
		    echo "beta=$beta" >> $dir/config.sh
		    echo "szbeam=1.2" >> $dir/config.sh
		    echo "sz_xi=$spt_xi" >> $dir/config.sh



		    if [ "$center" = "xrayNONE" ]; then
			echo "xraycentering=None" >> $dir/config.sh

		    elif [ "$center" = "xraymag" ]; then
			echo "xraycentering=magneticum" >> $dir/config.sh

		    elif [ "$center" = "szanalytic" ]; then
			echo "sztheoreticalcentering=analytic" >> $dir/config.sh

		    elif [ "$center" = "szmag" ]; then
			echo "sztheoreticalcentering=magneticum" >> $dir/config.sh

		    fi




		done
   		
	    done
	    	
	done
    
    }

done
	
