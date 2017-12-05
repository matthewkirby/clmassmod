#!/bin/bash -u

#profile fit range test

# Test run to see if new mass-concentration fitting will work.

rootdir=/project/kicp/avestruz/storage/mxxlsims/mxxlsnap
snap=54 # z=1 snapshot (54 is the z=0.25 snapshot)
rm ../run30mxxl$snap
# Want a radial range that encompasses the acs snapshot (3.2 arcmin on
# a side).  We want to exclude very center, but always want to go as
# far out as possible
# radius that matches ACS FoV
noise=n2_4

sampler=samplemcmc

# Note: core_none is perfect center, e.g. no offset
for center in core_none; do
    for r in rACS r10 ; do 
    for mc in cfree; do

	config=acs-${mc}-${r}-${center}-${noise}-${sampler}-aug2017

	echo $config >> ../run30mxxl$snap
		
	dir=$rootdir$snap/$config

	if [ ! -e $dir ]; then
	    mkdir -p $dir
	fi

	# lineargaussbins12 is how we're binning shear profile, n2_4 is noise, reassigndensitypicker

	# using samplemcmc.py because we are running a chain for both
	# mass and concentration.  Alternatively can use scanpdf.py if
	# we're only looking at the mass pdf

	# Need to combine composite with the fov picker: combinegalaxypickers.py, combining acs_singlepoint and n2_4
	
	cat overhead.py linearprior.py nozscaling.py fixedsourceredshift.py lineargaussbins12.py ${noise}.py acs_singlepoint.py nobinnoise.py mxxl.py ${mc}.py ${sampler}.py ${r}.py ${center}.py combinegalaxypickers.py > $dir/config.py

	echo "zsource=2.0" >> $dir/config.py

	done
    done
done
