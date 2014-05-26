#!/usr/bin/env python
#######################
# This file takes the MXXL binary file and converts it into an easy to use python format

# For now, I'm just going to use the existing file, and not include the info I don't have from Aaron's file
#######################

import sys, binaryutils, readtxtfile, cPickle
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

    redshift = readtxtfile.readtxtfile('mxxl_imperial/snap%d/redshift' % snapnum)[0,0]

    halocat = asciireader.read('mxxl_imperial/snap%d/halos_%d.txt' % (snapnum, snapnum))

    idconversion = readtxtfile.readtxtfile('mxxl_imperial/snap%d/idconversion' % snapnum)

    clusterinfo = {}

    try:
        for curindexs in [map(int, x) for x in idconversion]:

            cid, sid, projid = curindexs

            try:

                clusterinfo[cid] = dict(m500 = halocat['Mass500_crit'][sid]*1e10,
                                        m200 = halocat['Mass200'][sid]*1e10,
                                        concen = 1./halocat['EinastoParameter_1'][sid],
                                        redshift = redshift)

            except:

                clusterinfo[cid] = dict(m500 = 0.,
                                        m200 = halocat['Mass200'][sid]*1e10,
                                        concen = 0.,
                                        redshift = redshift)
                

    except ValueError:
        print curindexs

    with open(outfile, 'wb') as output:

        cPickle.dump(clusterinfo, output, -1)



        
