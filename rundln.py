#!/usr/bin/env python
###########################
# Run dln fit for one mass bin of one noise sim
###########################

import sys
import deconvolvedlognorm as dln
import nfwfit
import nfwutils
import numpy as np
import pymc

def defineMassEdges(simtype, delta):

    if simtype == 'bk11snap124' or simtype == 'bk11snap141':
        if delta == 200:
            massedges = np.array([4e14, 5e15])
        else:
            massedges = np.array([1.3e14, 5e15])

    elif simtype == 'mxxlsnap41':
        if delta == 200:
            massedges = np.logspace(np.log10(3.5e14), np.log10(1.4e15), 8)
        elif delta == 500:
            massedges = np.logspace(np.log10(3.2e14), np.log10(1.2e15), 5)
        elif delta == 2500:
            massedges = np.array([1.6e14, 4e14])

    elif simtype == 'mxxlsnap54':
        if delta == 200:
            massedges = np.logspace(np.log10(1.5e14), np.log10(4e15), 10)
        elif delta == 500:
            massedges = np.logspace(np.log10(7.8e14), np.log10(3e15), 7)
        elif delta == 2500:
            massedges = np.array([3.15e14, 1.2e15])


    return massedges

#########

    

def run(simtype, chaindir, outfile, delta, pdftype, massbin=0):

    config = nfwfit.readConfiguration('%s/config.sh' % chaindir)
    simreader = nfwfit.buildSimReader(config)
    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    massedges = defineMassEdges(simtype, delta)

    isPDF = False
    if pdftype == 'pdf':
        isPDF = True
        halos = dln.loadPosteriors(chaindir, simtype, simreader, delta, massedges, massbin,
                                   reader = dln.PDFReader)
    elif pdftype == 'mcmc':
        halos = dln.loadPosteriors(chaindir, simtype, simreader, delta, massedges, massbin,
                                   reader = dln.MCMCReader,
                                   cprior = 100.)
        
        

    if len(halos) < 10:
        sys.exit(0)

    success = False
    for i in range(20):

        try:
            
            if isPDF == True:
                parts = dln.buildPDFModel(halos)
            else:
                parts = dln.buildMCMCModel(halos)

            success=True
            break

        except (AssertionError, pymc.ZeroProbability) as e:
            
            continue

    assert(success is True)

    with open('%s.massrange' % outfile, 'w') as output:
        output.write('%f\n%f\n' % (massedges[massbin], massedges[massbin+1]))
    dln.sample(parts, outfile, 10000, singlecore=True)




if __name__ == '__main__':

    massbin = 0

    simtype=sys.argv[1]
    chaindir=sys.argv[2]
    outfile=sys.argv[3]
    delta=int(sys.argv[4])
    pdftype=sys.argv[5]
    if len(sys.argv) > 6:
        massbin=int(sys.argv[6])

    print 'Called with:', dict(simtype=simtype, chaindir=chaindir, 
                               outfile=outfile, delta=delta, 
                               pdftype = pdftype,
                               massbin=massbin)

    
    run(simtype, chaindir, outfile, delta, pdftype, massbin)
        
