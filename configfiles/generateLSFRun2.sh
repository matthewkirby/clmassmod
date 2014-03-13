#!/bin/bash

for i in 0 3 4; do

    for j in 0 1 2; do

	for r in 1 5 10; do

	    for mc in c4 cfree; do

		config=${mc}-r${r}-n${i}_${j}
		dir=~/nfs/bcc_clusters/recentered/$config
		mkdir $dir
		cat base50.sh bcc.sh ${mc}.sh r${r}.sh n${i}_${j}.sh > $dir/config.sh

		echo $config >> ../bccrun2

	    done

	done

    done

done
	