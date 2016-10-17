#!/bin/bash -u

snap=41
sampling=6_4


for binning in linearbins6 linearbins18 logbins6 logbins12 logbins18; do
    
    for r in 5 10; do
	
	for mc in c4 duffy; do
	
	    for coresize in 0 5 none; do
	    
		config=${mc}-r${r}-n${sampling}-core${coresize}-${binning}
		dir=../mxxl_imperial/mxxlsnap$snap/$config
		mkdir $dir
		cat maxlike.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
	    
		echo $config >> ../run9
	    
	    done

	done
   	
    done
    
done


	