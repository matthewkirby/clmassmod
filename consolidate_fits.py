#!/usr/bin/env python
############################

import glob, cPickle, sys, os, re
import numpy as np
import nfwutils, readMXXL

###########################

simdir=sys.argv[1]
outdir=sys.argv[2]

medians = []
sigma = []
actuals = []
redshifts = []

siminfo_filename = '{0}/siminfo.pkl'.format(simdir)

siminfo = cPickle.load(open(siminfo_filename, 'rb'))

redshift = siminfo['redshift']
massmapping = siminfo['massmapping']

haloid_pattern = re.compile('halo_cid(\d+)\.out')

class WeirdException(Exception): pass

for output in glob.glob('%s/*.out' % outdir):

    try:

        filebase = os.path.basename(output)
        
        match = haloid_pattern.match(filebase)

        haloid = int(match.group(1))

        
        actual_mass = massmapping[haloid][0]*1e10


        input = open(output)
        m200s, nfails = cPickle.load(input)
        input.close()

        if len(m200s) == 0:
            print 'All failed in {0}'.format(output)
            continue


        if len(m200s) == 1:
            median = m200s[0]
            sym_std = 0.

        else:
            sorted_m200s = np.sort(m200s)
            median = sorted_m200s[int(0.5*len(m200s))]
            r68low = median - sorted_m200s[int(0.16*len(m200s))]
            r68high = sorted_m200s[int(0.84*len(m200s))] - median
            sym_std = (r68high + r68low)/2.

        medians.append(median)
        sigma.append(sym_std)

        actuals.append(actual_mass)
        redshifts.append(redshift)

    except:

        print 'Skipping ', output


results = map(np.array, [medians, sigma, actuals, redshifts])

print 'Consolidated {0} fits'.format(len(medians))



cPickle.dump(results, open('%s/consolidated.pkl' % outdir, 'w'))

    
    
