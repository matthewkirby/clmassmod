#!/usr/bin/env python
#######################
# This file takes the MXXL binary file and converts it into an easy to use python format

# For now, I'm just going to use the existing file, and not include the info I don't have from Aaron's file
#######################

import sys, binaryutils, readtxtfile, cPickle

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

def readHalosFile(filename):

    data = readtxtfile.readtxtfile(filename)

    params = 'Mass200	Radius200	Position_0	Position_1	Position_2	VirialRatio	CenterOfMassOffset	SubstructureFraction	EinastoParameter_0	EinastoParameter_1	EinastoParameter_2	HalfMassRadius	HalfMassFormationRedshift	SubhaloFileNumber	SubhaloFileOffset	EinastoMeanRhoR2	MeanRhorh	EinastoFormationRedshift'.split()

    halocat = {}

    for i, name in enumerate(params):

        halocat[name] = data[:,i]

    return halocat
    

######################

if __name__ == '__main__':

    outfile = sys.argv[1]

    siminfo = cPickle.load(open('mxxl_imperial/snap41/siminfo.pkl'))
    massmapping = siminfo['massmapping']
    redshift = siminfo['redshift']

    halocat = readHalosFile('mxxl_imperial/snap41/halos_41.txt')

    idconversion = readtxtfile.readtxtfile('mxxl_imperial/snap41/idconversion')

    clusterinfo = {}

    try:
        for curindexs in [map(int, x) for x in idconversion]:

            cid, sid, projid = curindexs

            clusterinfo[cid] = dict(m500 = halocat['Mass500_crit'][sid]*1e10,
                                    m200 = halocat['Mass200'][sid]*1e10,
                                    concen = 1./halocat['EinastoParameter_1'][sid],
                                    redshift = redshift)

            assert(clusterinfo[cid]['m200'] == (massmapping[cid][0]*1e10))
    except ValueError:
        print curindexs

    with open(outfile, 'wb') as output:

        cPickle.dump(clusterinfo, output, -1)



        
