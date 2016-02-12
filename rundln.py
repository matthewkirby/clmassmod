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
import pymc
import cPickle
import simutils


def defineBins(masses, minclusters = 56, maxclusters = 600, binstep=0.2):

    if len(masses) < minclusters:
        return [np.min(masses), np.max(masses)]

    nmasses = len(masses)
    smass = np.sort(masses)[::-1]
    l10smass = np.log10(smass)
    indexes = np.arange(nmasses)
    edges = [np.log10(smass[0] + 1e12)]
    
    startindex = 0

    while(startindex < (nmasses-minclusters)):


        minbinsize = l10smass[startindex] - l10smass[startindex+minclusters]

        maxstopindex = np.min([startindex+maxclusters, nmasses-1])
        maxbinsize = l10smass[startindex] - l10smass[maxstopindex]
        ismax = False

        if minbinsize > binstep:
            print 'Min'
            stopval = l10smass[startindex+minclusters]
        elif maxbinsize < binstep:
            print 'Max'
            ismax = True
            stopval = l10smass[maxstopindex]
            if maxstopindex == (nmasses -1):
                break
        else:
            print 'Interval'
            stopval = l10smass[startindex]-binstep

        stopindex = indexes[l10smass < stopval][0]
        if (nmasses - stopindex) < minclusters or (ismax and (nmasses - stopindex < maxclusters)) :
            print 'End Adjust'
            equalstop = int((nmasses + startindex)/2.)
            stopval = l10smass[equalstop]
            stopindex = indexes[l10smass < stopval][0]

        edges.append(stopval)
        startindex = stopindex



    edges.append(l10smass[-1])

    return 10**np.array(edges[::-1])


class NoDataException(Exception): pass

def defineMassEdges(simtype, delta):

    if simtype == 'bk11snap124':
        if delta == 500:
            #3 bins
            massedges = np.array([  1.5e+14,   2.2e+14,   3.5e+14,6.5e+14])
        elif delta == 200:
            #3 bins
            massedges = np.array([  3.0e+14,   3.5e+14,   5.0e+14, 9.9e+14])
        else:
            raise NoDataException

    elif simtype == 'bk11snap141':
        if delta == 500:
            #4 bins
            massedges = np.array([2.0e+14, 2.3e+14, 2.9e+14, 4.4e+14, 8.8e+14])
        elif delta == 200:
            #2 bins
            massedges = np.array([4.5e+14,   6.7e+14,  1.6e+15])
        else:
            raise NoDataException

    elif simtype == 'mxxlsnap41':
        if delta == 200:
            #5 bins
            massedges = np.array([3.5e+14, 3.9e+14, 4.5e+14, 6.4e+14, 9.1e+14, 1.4e+15])
        elif delta == 500:
            #4 bins
            massedges = np.array([3.0e+14, 3.2e+14, 4.6e+14, 6.5e+14, 1.2e+15])
        elif delta == 2500:
            #3 bins
            massedges = np.array([1.6e+14, 1.8e+14, 2.4e+14, 3.9e+14])

    elif simtype == 'mxxlsnap54':
        if delta == 200:
            massedges = np.logspace(np.log10(1.5e14), np.log10(4e15), 10)
        elif delta == 500:
            massedges = np.logspace(np.log10(7.8e14), np.log10(3e15), 7)
        elif delta == 2500:
            massedges = np.array([3.15e14, 1.2e15])


    return massedges

#########

    

def run(simtype, chaindir, outfile, delta, pdftype, massbin=0, sigmapriorfile = None):

    config = simutils.readConfiguration('%s/config.py' % chaindir)
    simreader = config['simreader']
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

    if len(halos) > 500:
        print 'Down sampling'
        halos = random.sample(halos, 500)

    sigmapriors = None
    if sigmapriorfile is not None:
        print 'Using Sigma Priors!'
        sigmapriors = cPickle.load(open(sigmaprior, 'rb'))[massbin]



    success = False
    for i in range(20):

        try:
            
            if isPDF == True:
                parts = dln.buildPDFModel(halos, sigmapriors = sigmapriors)
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
    sigmaprior = None

    simtype=sys.argv[1]
    chaindir=sys.argv[2]
    outfile=sys.argv[3]
    delta=int(sys.argv[4])
    pdftype=sys.argv[5]
    if len(sys.argv) > 6:
        massbin=int(sys.argv[6])
    if len(sys.argv) > 7:
        sigmaprior = sys.argv[7]

    print 'Called with:', dict(simtype=simtype, chaindir=chaindir, 
                               outfile=outfile, delta=delta, 
                               pdftype = pdftype,
                               massbin=massbin,
                               sigmaprior = sigmaprior)

    
    run(simtype, chaindir, outfile, delta, pdftype, massbin, sigmaprior)
        
