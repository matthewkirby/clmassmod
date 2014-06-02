#!/bin/bash

for snap in 41 54; do

    for i in 0 3 4; do

	for r in 8 9 10; do

	    for mc in zhaomc duffy; do

		config=${mc}-r${r}-n${i}_0
		dir=../mxxl_imperial/snap$snap/$config
		mkdir $dir
		cat base.sh mxxl.sh ${mc}.sh r${r}.sh n${i}_0.sh > $dir/config.sh

		echo $config >> ../run5

	    done

	done

    done

done
	