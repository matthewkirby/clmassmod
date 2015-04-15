#!/bin/bash -u

## Explore the effect of having a density transition in the fit

snap=41
r=5
mc=c4

for targetz in 0.5 1.0; do
	    
    config=lineargaussbins-${mc}-r${r}-targetz${targetz}-splitdensity-sige28
    dir=../../mxxl_lensing/mxxlsnap$snap/$config
    mkdir $dir
    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh splitdensity.sh lineargaussbins12.sh > $dir/config.sh
    echo "targetz=${targetz}" >> $dir/config.sh
    echo "" >> $dir/config.sh
	    
    echo $config >> ../run12e
	    
done

for noise in 7_5 4_5; do

    config=lineargaussbins-${mc}-r${r}-targetz1.0-n${noise}
    dir=../../mxxl_lensing/mxxlsnap$snap/$config
    mkdir $dir
    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh lineargaussbins12.sh n${noise}.sh > $dir/config.sh
    echo "targetz=1.0" >> $dir/config.sh
    echo "" >> $dir/config.sh
    echo $config >> ../run12e

done
    

	
