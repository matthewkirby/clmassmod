#!/usr/bin/env python
#######################
# This file compiles the truth tables for BCC clusters
#######################

import sys, os, glob, cPickle, re
import astropy.io.fits as pyfits
import readtxtfile, ldac

#######################


simdir = sys.argv[1]
outfile = sys.argv[2]

clusterfiles = glob.glob('{0}/haloid*fit'.format(simdir))

clusterinfo = {}

idmatch = re.compile('haloid(\d+)_zLens.+')

for clusterfile in clusterfiles:

    base = os.path.basename(clusterfile)
    match = idmatch.match(base)
    if match is None:
        print 'Skipping {0}'.format(clusterfile)
        continue

    id = int(match.group(1))
    
    halocat = ldac.LDACCat(pyfits.open(clusterfile)[1])
    
    
    clusterinfo[id] = dict(m500 = halocat['M500C'][0],
                           m200 = halocat['M200C'][0],
                           concen = halocat['C200C'][0],
                           redshift = halocat['ZLENS'][0])





with open(outfile, 'wb') as output:

    cPickle.dump(clusterinfo, output, -1)





    

    


    
