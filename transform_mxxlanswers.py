#!/usr/bin/env python
#######################
# This file takes the MXXL binary file and converts it into an easy to use python format

# For now, I'm just going to use the existing file, and not include the info I don't have from Aaron's file
#######################

import sys, binaryutils, readtxtfile, cPickle, glob, os, re
import astropy.io.ascii as asciireader
import numpy as np

######################





def readMXXLBinary(filename):

    with open(filename, 'rb') as input:

        nelem = binaryutils.readVal(input, 'I')

        properties = {}

        varnames = 'm200 r200 x y z Vir doff fsub SubFile SubLoc r500_crit M500_crit'.split()

        for varname in varnames:
            properties[varname] = binaryutils.readArray(input, 'f', shape=(nelem,))



    return properties

########################


if __name__ == '__main__':

    snapnum = int(sys.argv[1])
    outfile = sys.argv[2]

    redshift = readtxtfile.readtxtfile('/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap%d/redshift' % snapnum)[0,0]

    halocat = readtxtfile.readtxtfile('/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap%d/halos_%d.txt' % (snapnum, snapnum))

    if snapnum == 41:
        properties = readMXXLBinary('/users/dapple/astro/mxxlsims/mxxl_imperial/mxxlsnap41/Doug_XXL_z1.bin')
        sort_props = np.argsort(properties['m200'])[::-1]

    clusterinfo = {}

    for halofile in glob.glob('/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap%d/*.convergence_map' % snapnum):

        halobase = os.path.basename(halofile)
        
        match = re.match('halo_%d_((\d+)_\d)\.convergence_map' % snapnum, 
                         halobase)
        
        myid = match.group(1)
        stefan_id = int(match.group(2))

        m500 = 0.
        if snapnum == 41:
            m500 = properties['M500_crit'][sort_props[stefan_id]]*1e10
       

        clusterinfo[myid] = dict(m500 = m500,
                                 m200 = halocat[stefan_id, 0]*1e10,
                                 concen = 0.,
                                 redshift = redshift)
                


    with open(outfile, 'wb') as output:

        cPickle.dump(clusterinfo, output, -1)



        
