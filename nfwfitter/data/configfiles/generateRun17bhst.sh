#!/bin/bash -u

snap=41
binning=hstnoisebins

cat ../shearprofiles/coresizeindex.list | { 

    while read cluster coresizeindex zcluster; do 
	
	for r in 5 16; do
	    
	    for mc in c4 duffy diemer15; do
		
		for center in xrayNONE xrayXVP szxvptcenter core${coresizeindex} xraylensingpeak szlensingpeak; do

		    config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    

		    if [ $center != "xrayXVP" ] && [ $center != "szxvptcenter" ]; then

			isdone=`grep $config /vol/euclid1/euclid1_2/dapple/rundlns/finished | wc -l`
			if [ "$isdone" -gt 0 ]; then
			    echo $config >> ../run17bhst.alreadyfinished
			    continue
			fi

			ishaloprocessed=`grep $config ../haloprocessed | wc -l`
			if [ "$ishaloprocessed" -gt 0 ]; then
			    echo $config >> ../run17bhst.needsdln
			    continue
			fi

		    fi

		    echo $config >> ../run17bhst.torun
		    

		    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config

		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		    
		    
		    cat scanpdf.sh mxxl.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
		    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
		    echo "targetz=$zcluster" >> $dir/config.sh
		    
		    #add center to config file
		    if [ "$center" == "szxvptcenter" ]; then
			echo "sztheoreticalcentering=xvp" >> $dir/config.sh
		    elif [ "$center" == "core${coresizeindex}" ]; then
			cat core_${coresizeindex}.sh >> $dir/config.sh
		    elif [ "$center" == "xrayXVP" ]; then
			echo "xraycentering=XVP" >> $dir/config.sh
		    elif [ "$center" == "xraylensingpeak" ]; then
			echo "xraycentering=lensingpeak" >> $dir/config.sh
		    elif [ "$center" == "szlensingpeak" ]; then
			echo "sztheoreticalcentering=lensingpeak" >> $dir/config.sh
		    fi
		    
		    
		done
   		
	    done
	    
	done
	
    done
    
}





	
