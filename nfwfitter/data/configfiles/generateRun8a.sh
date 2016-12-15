#!/bin/bash -u
############
# Use Gaussian Approx errors in the shear profile. Is this a problem with bootstrapping?
###################################

snap=41
binning=linearbins12

for sampling in 0_0 2_2 4_3 6_4; do
    
    for r in 5; do
	
	for mc in c4; do
	
	    for coresize in none; do
	    
		config=${mc}-r${r}-n${sampling}-core${coresize}-${binning}-gaussianshearerr
		dir=../mxxl_imperial/mxxlsnap$snap/$config
		mkdir $dir
		cat maxlike.sh mxxl_hstbeta.sh gaussianshearerrs.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
	    
		echo $config >> ../run8a
	    
	    done

	done
   	
    done
    
done


	