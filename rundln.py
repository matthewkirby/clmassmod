#!/usr/bin/env python
###########################
# Run dln fit for one mass bin of one noise sim
###########################

import sys
import deconvolvedlognorm as dln
import nfwfit
import nfwutils
import numpy as np

def run(simtype, chaindir, delta, massbin=0):

    config = nfwfit.readConfiguration('%s/config.sh' % chaindir)
    simreader = nfwfit.buildSimReader(config)
    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    if simtype == 'bk11snap124' or simtype == 'bk11snap141':
        if delta == 200:
            massedges = np.array([4e14, 5e15])
        else:
            massedges = np.array([1.3e14, 5e15])

    elif simtype == 'mxxlsnap41' or simtype == 'mxxlsnap54':
        massedges = np.logspace(np.log10(2e14), np.log10(3e15), 9)
        

    halos = dln.loadPDFs(chaindir, simtype, simreader, massedges, massbin)
    #halos = dln.loadPDFs(chaindir, simtype, simreader)

    if len(halos) < 5:
        sys.exit(0)

    parts = dln.buildPDFModel(halos)

    dln.sample(parts, '%s/dln_%d.%d' % (chaindir, massbin, delta), 10000, singlecore=True)




if __name__ == '__main__':

    massbin = 0

    simtype=sys.argv[1]
    chaindir=sys.argv[2]
    if len(sys.argv) == 5:
        massbin=int(sys.argv[3])
    delta=int(sys.argv[-1])

    print 'Called with:', simtype, chaindir, massbin, delta


    run(simtype, chaindir, delta, massbin)
