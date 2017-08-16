#!/bin/bash -u

#profile fit range test

# New science question to test with this run: What is the bias for
# single orbin HST-observations and how does that change with M-c
# relations, and X-ray or SZ centers?

rootdir=/project/kicp/avestruz/storage/mxxlsims/mxxlsnap
snap=41 # z=1 snapshot (54 is the z=0.25 snapshot)
rm ../run26mxxl$snap
# Want a radial range that encompasses the acs snapshot (3.2 arcmin on
# a side).  We want to exclude very center, but always want to go as
# far out as possible
r=r_run26
noise=n2_4
# Note: core_none is perfect center, e.g. no offset
for center in xraymagneticum szmagneticum_ignorecore core_none; do

    for mc in c3 c4 c5 diemer15; do

	config=acs-${mc}-${r}-${center}-${noise}-aug2017

	echo $config >> ../run26mxxl$snap
		
	dir=$rootdir$snap/$config

	if [ ! -e $dir ]; then
	    mkdir -p $dir
	fi

	# lineargaussbins12 is how we're binning shear profile, n2_4 is noise, reassigndensitypicker
	# scanpdf.py is because we're only looking at the mass pdf, but we would change this if we want to run a chain.

	# Need to combine composite with the fov picker: combinegalaxypickers.py, combining acs_singlepoint and n2_4
	
	cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py ${noise}.py acs_singlepoint.py nobinnoise.py mxxl.py ${mc}.py scanpdf.py ${r}.py ${center}.py combinegalaxypickers.py > $dir/config.py

	echo "zsource=2.0" >> $dir/config.py


    done
done
