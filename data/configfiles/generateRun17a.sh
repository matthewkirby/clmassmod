#!/bin/bash -u

binning=linearbins12
shapenoise=0.25

for snap in 124 141; do


    #MEGACAM
    cat megacam_siminput.reduced.list | { 

	while read cluster zcluster ndensity beta core coreindex; do 

	    for r in 9; do

		for mc in c4 duffy diemer15; do

		    for center in sztcenter corenone szxvptcenter core${coreindex}; do
	  	    
			    config=mega-${mc}-r${r}-sigma${shapenoise}-${center}-${cluster}
			    dir=../../bk11_lensing/snap$snap/intlength400/$config


#			    #link to euclid1_1
#
#			    noutputs=`ls $dir | wc -l`
#			    if [ $noutputs -gt 0 ]; then
#				rundlndir=/vol/euclid1/euclid1_2/dapple/rundlns/bk11snap${snap}/$config
#				noutputs=`ls $rundlndir | wc -l`
#				if [ $noutputs -gt 0 ]; then
#				    echo $config >> ../run17a.${snap}.alreadydone
#				else
#				    echo $config >> ../run17a.${snap}.maybedone
#				fi
#				continue
#			    fi
#
			    realdir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
			    mkdir -p $realdir

#			    if [ -e $dir ]; then
#				mv $dir ../../mxxl_lensing/mxxlsnap$snap/old_${config}
#			    fi
#			    ln -s $realdir $dir

			    dir=$realdir
#

			    cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
			    echo "targetz=$zcluster" >> $dir/config.sh
			    echo "nperarcmin=$ndensity" >> $dir/config.sh
			    echo "shapenoise=$shapenoise" >> $dir/config.sh
			    echo "beta=$beta" >> $dir/config.sh

			    #add center to config file
			    if [ "$center" == "szxvptcenter" ]; then
				echo "sztheoreticalcentering=xvp" >> $dir/config.sh
			    elif [ "$center" == "core${coreindex}" ]; then
				cat core_${coreindex}.sh >> $dir/config.sh
			    elif [ $"center" == "sztcenter" ]; then
				echo "sztheoreticalcentering=True" >> $dir/config.sh
			    fi
			    
#			    echo $config >> ../run17a.${snap}
	       
			
		    done
   		    
		done
		
	    done
	    
	done
    
    }

done

#snap=124
#binning=hstnoisebins
#
#
#cat ../shearprofiles/coresizeindex.list | { 
#
#    while read cluster coresizeindex zcluster; do 
#	
#	for r in 5 16; do
#	    
#	    for mc in c4 duffy diemer15 bhattacharya13; do
#		
#		for center in xrayNONE xrayXVP szxvptcenter core${coresizeindex} xraylensingpeak szlensingpeak; do
#		    config=hstnoisebins-${mc}-r${r}-${center}-${cluster}	  	    
#		    dir=/vol/euclid1/euclid1_raid1/dapple/bk11_lensing/snap$snap/intlength400/$config
#		    realdir=/vol/euclid1/euclid1_1/dapple/bk11_lensing/snap$snap/intlength400/$config
#		    
#		    
##		    #link to euclid1_1
##		    if [ -e $dir ]; then
##			noutputs=`ls $dir | wc -l`
##			if [ $noutputs -gt 1 ]; then
##			    rundlndir=/vol/euclid1/euclid1_2/dapple/rundlns/bk11snap${snap}/$config
##			    noutputs=`ls $rundlndir | wc -l`
##			    if [ -e $rundlndir] && [ $noutputs -gt 0 ]; then
##				echo $config >> ../run17a.${snap}.alreadydone
##			    else
##				echo $config >> ../run17a.${snap}.maybedone
##			    fi
##			    continue
##
##			fi
##			
##		    fi
##
##
##		    
##
##		    mkdir $realdir
##		    
##		    ln -s $realdir $dir
##		    
#		    
#		    cat scanpdf.sh bk11.sh ${mc}.sh r${r}.sh ${binning}.sh > $dir/config.sh
#		    echo "profilefile=/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/${cluster}.szcenter.profile" >> $dir/config.sh
#		    echo "targetz=$zcluster" >> $dir/config.sh
#		    
#		    #add center to config file
#		    if [ "$center" == "szxvptcenter" ]; then
#			echo "sztheoreticalcentering=xvp" >> $dir/config.sh
#		    elif [ "$center" == "core${coresizeindex}" ]; then
#			cat core_${coresizeindex}.sh >> $dir/config.sh
#		    elif [ "$center" == "xrayXVP" ]; then
#			echo "xraycentering=XVP" >> $dir/config.sh
#		    elif [ "$center" == "xraylensingpeak" ]; then
#			echo "xraycentering=lensingpeak" >> $dir/config.sh
#		    elif [ "$center" == "szlensingpeak" ]; then
#			echo "sztheoreticalcentering=lensingpeak" >> $dir/config.sh
#		    fi
#		    
##		    echo $config >> ../run17a.${snap}
#		    
#		    
#		done
#   		
#	    done
#	    
#	done
#	
#    done
#    
#}
#




	
