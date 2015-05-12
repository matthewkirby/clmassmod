#!/usr/bin/env python
#######################
# This file takes the MXXL binary file and converts it into an easy to use python format

# For now, I'm just going to use the existing file, and not include the info I don't have from Aaron's file
#######################

import sys, readtxtfile, cPickle, glob, os, re
import astropy.io.ascii as asciireader
import numpy as np
import readMXXLProfile as rmp

######################




snapnum = int(sys.argv[1])
outfile = sys.argv[2]

redshift = readtxtfile.readtxtfile('/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap%d/redshift' % snapnum)[0,0]

clusterinfo = {}

for halofile in glob.glob('/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap%d/*.convergence_map' % snapnum):

    halobase = os.path.basename(halofile)

    match = re.match('halo_%d_((\d+)_\d)\.convergence_map' % snapnum, 
                     halobase)

    myid = match.group(1)
    stefan_id = int(match.group(2))

    profile = rmp.MXXLProfile('/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap%d/halo_%d_%d.radial_profile_3D.txt' % (snapnum, snapnum, stefan_id))

    clusterinfo[myid] = dict(m500 = profile.overdensityMass(500.),
                             m200 = profile.overdensityMass(200.),
                             m1p5 = profile.massEnclosed(1.5),
                             m2500 = profile.overdensityMass(2500.),
                             concen = 0.,
                             redshift = redshift)



with open(outfile, 'wb') as output:

    cPickle.dump(clusterinfo, output, -1)



        
