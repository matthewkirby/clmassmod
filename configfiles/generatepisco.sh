#!/bin/bash



for setup in subaru_z4 pisco_z4_3 pisco_z4_4 pisco_z4_4square; do

	for mc in c4 cfree; do

	    config=${mc}-sour_${setup}
	    mkdir ../mxxl_imperial/snap41/$config
	    cat base.sh ${mc}.sh r10.sh ${setup}.sh > ../mxxl_imperial/snap41/$config/config.sh

	    echo $config >> ../sour

	done



done

