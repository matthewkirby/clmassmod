#!/bin/bash -u


binning=linearbins12
shapenoise=0.25
r=9

for snap in 124 141; do

    cat megacam_siminput.reduced.list | awk '($1 !~ "#"){print}' | { 

	while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
	
	    for mc in duffy diemer15; do

		for center in szanalytic szxvpbcg; do
		
		    config=mega-${mc}-r${r}-sigma${shapenoise}-${center}-${cluster}
		
		    echo $config >> ../run17ebk11$snap
		
		    dir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		    
		    
		    cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh core_${coresizeindex}.sh > $dir/config.sh
		    
		    echo "targetz=$zcluster" >> $dir/config.sh
		    echo "nperarcmin=$ndensity" >> $dir/config.sh
		    echo "shapenoise=$shapenoise" >> $dir/config.sh
		    echo "beta=$beta" >> $dir/config.sh

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

    }
    
done





	
