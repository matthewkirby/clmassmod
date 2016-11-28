#!/bin/bash -u

#profile fit range test


center=xrayNONE
mc=diemer15
snap=54
for r in 3 5 6 7 9 10 19 20; do

    config=general-${mc}-r${r}-${center}-n2_4-nov2016
		
    echo $config >> ../run24mxxl$snap
		
    dir=/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap$snap/$config
		
    if [ ! -e $dir ]; then
	mkdir $dir
    fi

		 
    cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py n2_4.py reassigndensitypicker.py nobinnoise.py mxxl.py ${mc}.py scanpdf.py r${r}.py core_none.py > $dir/config.py

    echo "zsource=2.0" >> $dir/config.py


done
	   	

#cosmo & redshift, mc, center

for sim in mxxl bk11; do

    if [ $sim == 'mxxl' ]; then
	snaps="41 54"
    elif [ $sim == 'bk11' ]; then
	snaps="124 141"
    fi

    for snap in $snaps; do

	echo "Current: $sim $snap"

	if [ $sim == 'mxxl' ]; then
	    subdir=mxxl_lensing/mxxlsnap$snap
	elif [ $sim == 'bk11' ]; then
	    subdir=bk11_lensing/snap$snap/intlength400
	fi
	

	for mc in diemer15 c4 duffy; do

	    for r in 6 10 19; do

		for center in xrayNONE szmag xraymag; do

		    config=general-${mc}-r${r}-${center}-n2_4-nov2016
		
		    echo $config >> ../run24${sim}${snap}
		
		    dir=/vol/euclid1/euclid1_1/dapple/$subdir/$config
		
		    if [ ! -e $dir ]; then
			mkdir $dir
		    fi


		    cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py n2_4.py reassigndensitypicker.py nobinnoise.py ${sim}.py ${mc}.py scanpdf.py r${r}.py  > $dir/config.py

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
    done
done
