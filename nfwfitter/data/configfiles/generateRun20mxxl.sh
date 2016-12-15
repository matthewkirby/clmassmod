#!/bin/bash -u


binning=linearbins12
noise=n0_0

for snap in 41 54; do
	
    for mc in c4 diemer15; do

	for r in 6 7 9 10 17 18; do

	    for center in xrayNONE xrayWTG; do
		
		config=general-${mc}-r${r}-${noise}-${center}
		
		echo $config >> ../run18bmxxl$snap
		
		dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		if [ ! -e $dir ]; then
		    mkdir $dir
		fi
		    
		    
		    
		cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh ${noise}.sh > $dir/config.sh
		    

		if [ "$center" == "xrayWTG" ]; then
			echo "xraycentering=WTG" >> $dir/config.sh
		fi
			
	     		    
   		
	    done
	    
	done

    done
    
done





	
