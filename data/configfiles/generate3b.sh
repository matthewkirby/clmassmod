#!/bin/bash

for i in 0 1 2 3 4; do

    for hst in deepmosaic cvzsquare; do

	for mc in c4 cfree duffy; do

	    config=${mc}-hst_${hst}_${i}
	    mkdir ../mxxl_imperial/snap41/$config
	    cat base50.sh ${mc}.sh r6.sh hst_${hst}_${i}.sh > ../mxxl_imperial/snap41/$config/config.sh

	    echo $config >> ../mxxlrun3b

	done

    done

done

