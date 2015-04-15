#!/bin/bash -u

## Explore the effect of having a density transition in the fit

snap=41
r=5
mc=c4

for targetz in 0.5 1.0; do

    for noise in 7_4 4_4 5_4; do

	config=lineargaussbins-${mc}-r${r}-targetz${targetz}-n${noise}
	dir=../../mxxl_lensing/mxxlsnap$snap/$config
	mkdir $dir
	cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh lineargaussbins12.sh n${noise}.sh > $dir/config.sh
	echo "targetz=${targetz}" >> $dir/config.sh
	echo "" >> $dir/config.sh
	echo $config >> ../run12a

	    
    done

done


	
