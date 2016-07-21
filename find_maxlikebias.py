#!/usr/bin/env python
###########################
# Run dln fit for one mass bin of one noise sim
###########################

import sys
import random
import deconvolvedlognorm as dln
import nfwfit
import nfwutils
import numpy as np
import cPickle
import simutils
import rundln


########    

def run(simtype, chaindir, outfile, delta, massbin=0):

    config = simutils.readConfiguration('%s/config.py' % chaindir)
    simreader = config['simreader']
    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())
    model = config['model']

    selectors = rundln.defineMassEdges(simtype, delta)
    selector = selectors[massbin]


    isPDF = True
    halos = dln.loadPosteriors(chaindir, simtype, simreader, delta, selector,
                               reader = dln.PDFReader, model = model)
        

    if len(halos) < 10:
        sys.exit(0)


    maxlikes = np.zeros(len(halos))
    truths = np.zeros_like(maxlikes)
    for i in range(len(halos)):
        curhalo = halos[i]
        
        maxlikes[i] = curhalo['masses'][curhalo['pdf'] == np.max(curhalo['pdf'])][0]
        truths[i] = curhalo['true_mass']

    ratios = maxlikes/truths
    
    with open('{}.maxlike'.format(outfile), 'w') as output:
        output.write('mean stddev meanerr\n')
        output.write('{} {} {}\n'.format(np.mean(ratios), np.std(ratios, ddof=1), np.std(ratios,ddof=1)/np.sqrt(len(halos))))

if __name__ == '__main__':

    massbin = 0

    simtype=sys.argv[1]
    chaindir=sys.argv[2]
    outfile=sys.argv[3]
    delta=int(sys.argv[4])
    if len(sys.argv) > 5:
        massbin=int(sys.argv[5])

    print 'Called with:', dict(simtype=simtype, chaindir=chaindir, 
                               outfile=outfile, delta=delta, 
                               massbin=massbin)

    
    run(simtype, chaindir, outfile, delta, massbin)
        
