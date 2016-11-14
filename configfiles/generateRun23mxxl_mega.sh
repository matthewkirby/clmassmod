binning=lineargaussbins12
shapenoise=0.25
r=9

for snap in 41 54; do

    cat megacam_siminput.list | awk '($1 !~ "#"){print}' | { 

	while read cluster zcluster ndensity beta core coresizeindex spt_xi; do 
		 
		
		for mc in c4 diemer15; do

		    config=mega-${mc}-r${r}-${center}-${cluster}-feb2016
		
		    echo $config >> ../run22bmxxl$snap
		
		    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi
		    
		   
		 
		    cat megacam.py mxxl.py ${mc}.py scanpdf.py r${r}.py > $dir/config.py


		    echo "targetz=$zcluster" >> $dir/config.py
		    echo "nperarcmin=$ndensity" >> $dir/config.py
		    echo "shapenoise=$shapenoise" >> $dir/config.py
		    echo "beta=$beta" >> $dir/config.py


		    if [ "$center" = "xrayNONE" ]; then
			cat core_none.py >> $dir/config.py

		    elif [ "$center" = "xraymag" ]; then
			cat xraymagneticum.py >> $dir/config.py

		    elif [ "$center" = "szanalytic" ]; then
			cat szanalytic.py >> $dir/config.py
			cat core_${coresizeindex}.py >> $dir/config.py
			echo "szbeam=1.2" >> $dir/config.py
			echo "sz_xi=$spt_xi" >> $dir/config.py

		    elif [ "$center" = "szmag" ]; then
			cat szmagneticum.py >> $dir/config.py
			cat core_${coresizeindex}.py >> $dir/config.py

		    fi


		done
   	       
	    	
	done
    
    }

done
