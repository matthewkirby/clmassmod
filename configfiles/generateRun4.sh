#!/bin/bash

## Set up test of different binning strategies

for binning in equalbins50 equalbins200 logbins6 logbins12 logbins18 linearbins6 linearbins12 linearbins18; do

    for mc in c4 cfree; do

	for density in n0_0 n3_0; do
	    
	    for rrange in r11 r5 r6 r10; do

		config=${mc}-${rrange}-${density}-${binning}
		
		mkdir ../mxxl_imperial/snap41/$config
		cat mxxl.sh $binning.sh ${mc}.sh $rrange.sh $density.sh  > ../mxxl_imperial/snap41/$config/config.sh

	    echo $config >> ../run4

	    done

	done

    done

done

