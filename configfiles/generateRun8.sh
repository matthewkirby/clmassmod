#!/bin/bash -u

snap=41
binning=linearbins12

for sampling in 0_0 2_2 4_3 6_4; do
    
    for r in 5 8; do
	
	for mc in c4 duffy; do
	
	    for coresize in 0 5 none; do
	    
		config=${mc}-r${r}-n${sampling}-core${coresize}-${binning}
		dir=../mxxl_imperial/mxxlsnap$snap/$config
		mkdir $dir
		cat base.sh maxlike.sh mxxl.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
	    
		echo $config >> ../run8
	    
	    done

	done
   	
    done
    
done


	