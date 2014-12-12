#!/usr/bin/env python
#######################
# This file takes the MXXL binary file and converts it into an easy to use python format

# For now, I'm just going to use the existing file, and not include the info I don't have from Aaron's file
#######################

import sys, binaryutils, readtxtfile, cPickle, glob, os, re
import astropy.io.ascii as asciireader

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

    redshift = readtxtfile.readtxtfile('/users/dapple/astro/mxxlsims/mxxl_imperial/mxxlsnap%d/redshift' % snapnum)[0,0]

    halocat = readtxtfile.readtxtfile('/users/dapple/astro/mxxlsims/mxxl_imperial/mxxlsnap%d/halos_%d.txt' % (snapnum, snapnum))

    clusterinfo = {}

    for halofile in glob.glob('/users/dapple/astro/mxxlsims/mxxl_imperial/mxxlsnap%d/*.convergence_map' % snapnum):

        halobase = os.path.basename(halofile)
        
        match = re.match('halo_((\d+)_\d)\.convergence_map', halobase)
        
        myid = match.group(1)
        stefan_id = int(match.group(2))
        

        clusterinfo[myid] = dict(m500 = 0.,
                                 m200 = halocat[stefan_id, 0]*1e10,
                                 concen = 0.,
                                 redshift = redshift)
                


    with open(outfile, 'wb') as output:

        cPickle.dump(clusterinfo, output, -1)



        
