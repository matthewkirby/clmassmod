#!/bin/bash -u

snap=41
binning=lineargaussbins12

for sampling in 0_0 2_2 4_3 6_4; do

    r=5
    mc=c4
    coresize=none
	

	    
    config=${mc}-r${r}-n${sampling}-core${coresize}-${binning}
    dir=../mxxl_imperial/mxxlsnap$snap/$config
    mkdir $dir
    cat maxlike.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
	    
    echo $config >> ../run10
	    


done
   	


	