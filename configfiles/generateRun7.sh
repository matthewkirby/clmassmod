#!/bin/bash

snap=41

for sampling in 0_0 3_2 4_3; do
    
    for r in 5 6 8 9; do
	
	mc=c4
	
	for coresize in 0 5 11 none; do
	    
	    config=${mc}-r${r}-n${sampling}_core${coresize}
	    dir=../mxxl_imperial/mxxlsnap$snap/$config
	    mkdir $dir
	    cat base.sh maxlike.sh mxxl.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh > $dir/config.sh
	    
	    echo $config >> ../run7
	    
	done
   	
    done
    
done


	