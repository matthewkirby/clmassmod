#!/bin/bash -u

snap=41
sampling=0_0
binning=hstnoisebins
r=5

cat ../shearprofiles/coresizeindex.list | { while read cluster coresizeindex; do 

	for mc in c4 duffy; do

	    for coresize in none ${coresizeindex}; do
	    
		config=hstnoisebins-${mc}-r${r}-core${coresize}-${cluster}
		dir=../mxxl_imperial/mxxlsnap$snap/$config
		mkdir $dir
		cat scanpdf.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
		echo "profilefile=/vol/braid1/vol1/dapple/mxxl/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
	    
		echo $config >> ../run11
	    
	    done

	done
   	
    done
    
}


	
