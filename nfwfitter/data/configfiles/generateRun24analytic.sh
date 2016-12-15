#!/bin/bash -u

#profile fit range test

targetdir=$1

center=xrayNONE
mc=diemer15
snap=54
for r in 3 5 7 9 20; do

    config=general-${mc}-r${r}-${center}-n2_4-nov2016
		
    echo $config >> ../run24analytic
		
    dir=$targetdir/$config
		
    if [ ! -e $dir ]; then
	mkdir $dir
    fi

		 
    cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py n2_4.py reassigndensitypicker.py nobinnoise.py analytic.py ${mc}.py scanpdf.py r${r}.py core_none.py > $dir/config.py

    echo "zsource=2.0" >> $dir/config.py


done
	   	

#cosmo & redshift, mc, center



for mc in diemer15 c4 duffy; do

    for r in 6 10 19; do

	for center in xrayNONE szmag xraymag; do

	    config=general-${mc}-r${r}-${center}-n2_4-nov2016
		
	    echo $config >> ../run24analytic
		
	    dir=$targetdir/$config
		
	    if [ ! -e $dir ]; then
		mkdir $dir
	    fi


	    cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py n2_4.py reassigndensitypicker.py nobinnoise.py analytic.py ${mc}.py scanpdf.py r${r}.py  > $dir/config.py

	    echo "zsource=2.0" >> $dir/config.py
		    
	    if [ "$center" = "xrayNONE" ]; then
		cat core_none.py >> $dir/config.py

	    elif [ "$center" = "xraymag" ]; then
		cat xraymagneticum.py >> $dir/config.py

	    elif [ "$center" = "szmag" ]; then
		cat szmagneticum_ignorecore.py >> $dir/config.py
	    fi

 
	done
    done
done


