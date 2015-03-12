#!/usr/bin/env python
###########################
# Run dln fit for one mass bin of one noise sim
###########################

import sys
import deconvolvedlognorm as dln
import nfwfit
import numpy as np

simtype=sys.argv[1]
chaindir=sys.argv[2]
massbin=int(sys.argv[3])
#massbin='all'

config = nfwfit.readConfiguration('%s/config.sh' % chaindir)
simreader = nfwfit.buildSimReader(config)
massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)

halos = dln.loadPDFs(chaindir, simtype, simreader, massedges, massbin)
#halos = dln.loadPDFs(chaindir, simtype, simreader)

if len(halos) < 5:
    sys.exit(0)

parts = dln.buildPDFModel(halos)

dln.sample(parts, '%s/dln_%d' % (chaindir, massbin), 10000, singlecore=True)

