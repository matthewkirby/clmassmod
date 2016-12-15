#!/bin/bash -u


binning=linearbins12
noise=n0_0

for snap in 124 141; do
	
    for mc in c4 diemer15; do

	for r in 6 7 9 10 17 18; do

	    for center in xrayNONE xrayWTG; do
		
		config=general-${mc}-r${r}-${noise}-${center}
		
		echo $config >> ../run18bk11$snap
		
		dir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
		
		if [ ! -e $dir ]; then
		    mkdir $dir
		fi
		    
		    
		    
		cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh ${noise}.sh > $dir/config.sh
		    

		if [ "$center" == "xrayWTG" ]; then
			echo "xraycentering=WTG" >> $dir/config.sh
		fi
			
	     		    
   		
	    done
	    
	done

    done
    
done





	
