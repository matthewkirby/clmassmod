#!/bin/bash

snap=141

for sampling in 0_0 3_2 4_3; do
    
    for r in 6 7 9 10; do
	
	mc=c4
	
	for coresize in 0 5 11 none; do

	    config=${mc}-r${r}-n${sampling}_core${coresize}
	    dir=../../snap$snap/intlength400/$config
	    mkdir $dir
	    cat base.sh maxlike.sh bk11.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh > $dir/config.sh
	    
	    echo $config >> ../run7b

	    
	    
	done
   	
    done
    
done


