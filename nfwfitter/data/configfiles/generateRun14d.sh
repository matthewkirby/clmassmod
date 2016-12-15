#!/bin/bash -u

binning=linearbins12

rm ../run14b.*

for snap in 124 141; do

    cat megacam_siminput.list | { 

	while read cluster zcluster ndensity beta core coreindex; do 

	    for mc in c4 duffy; do
	    
		for r in 6 9; do
		    
		    for shapenoise in 0.25; do
	    
			    config=mega-${mc}-r${r}-sigma${shapenoise}-sztcenter-${cluster}
			    dir=../../bk11_lensing/snap$snap/intlength400/$config
			    mkdir $dir
			    cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
			    echo "targetz=$zcluster" >> $dir/config.sh
			    echo "nperarcmin=$ndensity" >> $dir/config.sh
			    echo "shapenoise=$shapenoise" >> $dir/config.sh
			    echo "beta=$beta" >> $dir/config.sh
			    echo "sztheoreticalcentering=True" >> $dir/config.sh
			    
			    
			    echo $config >> ../run14b.${snap}
	       
			
		    done
   		    
		done
		
	    done
	    
	done
    
    }

done


	
