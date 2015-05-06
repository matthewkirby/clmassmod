#!/usr/bin/env python
###########################
# Run dln fit for one mass bin of one noise sim
###########################

import sys
import deconvolvedlognorm as dln
import nfwfit
import nfwutils
import numpy as np

def run(simtype, chaindir, delta):

    config = nfwfit.readConfiguration('%s/config.sh' % chaindir)
    simreader = nfwfit.buildSimReader(config)
    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    if delta == 200:
        massedges = np.array([4e14, 5e15])
        massbin = 0
    else:
        massedges = np.array([1.3e14, 5e15])
        massbin=0
    #massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)

    halos = dln.loadPDFs(chaindir, simtype, simreader, massedges, massbin)
    #halos = dln.loadPDFs(chaindir, simtype, simreader)

    if len(halos) < 5:
        sys.exit(0)

    parts = dln.buildPDFModel(halos)

    dln.sample(parts, '%s/dln_%d.%d' % (chaindir, massbin, delta), 10000, singlecore=True)




if __name__ == '__main__':

    simtype=sys.argv[1]
    chaindir=sys.argv[2]
    delta=int(sys.argv[3])



    #massbin=int(sys.argv[3])
    #massbin='all'

    run(simtype, chaindir, delta)
