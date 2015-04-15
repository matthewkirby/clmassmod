#!/bin/bash -u

snap=54
binning=linearbins12

cat megacam_siminput.list | { 

    while read cluster zcluster ndensity beta core coreindex; do 

	for mc in c4 duffy; do

	    for r in 6 9; do

		for shapenoise in 0.25; do

		    for coresize in $coreindex none; do
	    
			config=mega-${mc}-r${r}-sigma${shapenoise}-core${coresize}-${cluster}
			dir=../../mxxl_lensing/mxxlsnap$snap/$config
			mkdir $dir
			cat scanpdf.sh mxxl_hstbeta.sh ${mc}.sh r${r}.sh core_${coresize}.sh ${binning}.sh > $dir/config.sh
			echo "targetz=$zcluster" >> $dir/config.sh
			echo "nperarcmin=$ndensity" >> $dir/config.sh
			echo "shapenoise=$shapenoise" >> $dir/config.sh
			echo "beta=$beta" >> $dir/config.sh
			
			
			echo $config >> ../run14
		    
		    done

		done
   	
	    done

	done

    done
    
}


	
