#!/usr/bin/env python
############################

import glob, h5py, h5pyutils, cPickle, sys, os
import numpy as np, nfwutils

###########################


masspklfile=sys.argv[1]
outdir=sys.argv[2]

medians = []
sigma = []
actuals = []
redshift = []


for output in glob.glob('%s/*.out' % outdir):

    try:

        adir, filebase = os.path.split(output)
        fileroot, fileext = os.path.splitext(filebase)

        simfile = '%s/%s.hdf5' % (simdir, fileroot)

        rawsim = h5py.File(simfile)
        shearcat = h5pyutils.getShearCat(rawsim)
        ra, dec, z, mass = h5pyutils.getClusterProperties(rawsim)
        rawsim.close()

        input = open(output)
        m200s, nfails = cPickle.load(input)
        input.close()

        sorted_m200s = np.sort(m200s)


        median = sorted_m200s[int(0.5*len(m200s))]
        r68low = median - sorted_m200s[int(0.16*len(m200s))]
        r68high = sorted_m200s[int(0.84*len(m200s))] - median
        sym_std = (r68high + r68low)/2.

        medians.append(median)
        sigma.append(sym_std)

        actuals.append(mass)
        redshift.append(z)

    except:

        print 'Skipping ', output


results = map(np.array, [medians, sigma, actuals, redshift])



cPickle.dump(results, open('%s/consolidated.pkl' % outdir, 'w'))

    
    
