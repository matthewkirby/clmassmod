#!/bin/bash -u

mc=cfree

for snap in 41 54; do
	
    for noise in n0_0 n3_4 n5_5; do

	if [ "${noise}" = "n0_0" ]; then
	    binning=linearbins12
	else
	    binning=lineargaussbins12
	fi

	for r in 5 9 10; do

	    for center in xrayNONE xrayWTG; do
		
		config=general-${mc}-r${r}-${noise}-${center}
		
		echo $config >> ../run19mxxl$snap
		
		dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		if [ ! -e $dir ]; then
		    mkdir $dir
		fi
		    
		    
		    
		cat mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh ${noise}.sh linearprior.sh > $dir/config.sh
		    

		if [ "$center" == "xrayWTG" ]; then
			echo "xraycentering=WTG" >> $dir/config.sh
		fi
			
	     		    
   		
	    done
	    
	done

    done
    
done





	
