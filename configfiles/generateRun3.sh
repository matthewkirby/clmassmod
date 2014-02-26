#!/bin/bash

for i in 0 1 2 3 4; do

    for hst in mosaic single square; do

	for mc in c4 cfree duffy; do

	    config=${mc}-hst_${hst}_${i}
	    mkdir ../mxxl_imperial/snap41/$config
	    cat base50.sh ${mc}.sh r6.sh hst_${hst}_${i}.sh > ../mxxl_imperial/snap41/$config/config.sh

	    echo $config >> ../run2

	done

    done

done


for i in 0 1 2 3 4 5; do

    for j in 0 1 2 3; do

	for r in 1 5 10; do

	    for mc in c4 cfree; do

		config=${mc}-r${r}-n${i}_${j}
		dir=../mxxl_imperial/snap41/$config
		mkdir $dir
		cat base50.sh mxxl.sh ${mc}.sh r${r}.sh n${i}_${j}.sh > $dir/config.sh

		echo $config >> ../run2

	    done

	done

    done

done
	