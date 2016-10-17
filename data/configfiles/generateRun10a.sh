#!/bin/bash -u

snap=41
binning=lineargaussbins12

for sampling in 2_2 4_3 6_4; do

    r=5
    mc=c4
    coresize=none
	

	    
    config=mcmc_linear-${mc}-r${r}-n${sampling}-core${coresize}-${binning}
    dir=../mxxl_imperial/mxxlsnap$snap/$config
    mkdir $dir
    cat mxxl_hstbeta.sh linearprior.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
	    
    echo $config >> ../run10a
	    
done

sampling=0_0
r=5
mc=c4
coresize=none
binning=linearbins12

	    
config=mcmc_linear-${mc}-r${r}-n${sampling}-core${coresize}-${binning}
dir=../mxxl_imperial/mxxlsnap$snap/$config
mkdir $dir
cat mxxl_hstbeta.sh linearprior.sh ${mc}.sh r${r}.sh n${sampling}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh

echo $config >> ../run10a




   	


	