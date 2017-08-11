#!/bin/bash -u

#profile fit range test


center=xrayNONE
mc=diemer15
snap=54
for r in 3 5 6 7 9 10 19 20; do

    for mc in c3 c4 c5; do

	config=general-${mc}-r${r}-${center}-n2_4-nov2016
		
	echo $config >> ../run25mxxl$snap
		
	dir=/home/dapple/storage/mxxlsims/mxxlsnap$snap/$config
		
	if [ ! -e $dir ]; then
	    mkdir $dir
	fi

		 
    cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py n2_4.py reassigndensitypicker.py nobinnoise.py mxxl.py ${mc}.py scanpdf.py r${r}.py core_none.py > $dir/config.py

    echo "zsource=2.0" >> $dir/config.py


    done
done

    
center=xraymag
mc=c4
snap=54
for r in 3 5 6 7 9 10 19 20; do

    config=general-${mc}-r${r}-${center}-n2_4-nov2016
		
    echo $config >> ../run25mxxl$snap
		
    dir=/home/dapple/storage/mxxlsims/mxxlsnap$snap/$config
		
    if [ ! -e $dir ]; then
	mkdir $dir
    fi

		 
    cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py n2_4.py reassigndensitypicker.py nobinnoise.py mxxl.py ${mc}.py scanpdf.py r${r}.py xraymagneticum.py > $dir/config.py

    echo "zsource=2.0" >> $dir/config.py


done
