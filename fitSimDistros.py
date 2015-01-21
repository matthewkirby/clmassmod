#!/usr/bin/env python
################################
# Take a sim fit consolidation file, and precompute lognormal intrinsic scatter models for mass bins
################################

import sys, os
import deconvolvedlognorm as dln
import pymc
import cPickle


def runFit(consol, inbin, outputfile, nsamples = 3000):

    parts = None
    for i in range(10):
        try:
            parts = dln.buildModel(measuredmass[inbin], measuredmasserr[inbin], truemass[inbin])
            break
        except pymc.ZeroProbability:
            continue
    if parts is None:
        raise pymc.ZeroProbability
    dln.sample(parts, outputfile, nsamples)



if __name__ == '__main__':

    consolfile = sys.argv[1]
    massbin = int(sys.argv[2])
    base, ext = os.path.splitext(consolfile)
    outputfile = '%s.lognorm.%d' % (base, massbin)

    with open(consolfile, 'rb') as input:
        data = cPickle.load(input)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)

    inbin = np.logical_and(data['true_m200s'] >= massedges[massbin],
                           data['true_m200s'] < massedges[massbin+1])

    if len(data['true_m200s'][inbin]) > 25:
        runFit(data, inbin, outputfile)

