#!/usr/bin/env python
####################
# Python MCMC program that uses pre-sampled M-C probabilities for individual clusters 
#  to measure the bias and M-C relation for lensing
#
#  Main approach is to importance sample the integrals needed, and to use the multiprocessing module for speed
####################

import numpy as np
import pymc


####################

def loadClusterData(simdir, chaindir):
    # loads M-C Chains for individual clusters
    pass

####################

def massbinModel(parts, clusters, massbinedges = np.logpsace(np.log10(1e14), np.log10(5e15), 10)):
    #constant spaced log mass bins

    nbins = len(massbinedges) - 1

    parts['logmassbiasratios'] = np.empty(nbins, dtype=object)
    parts['c200s'] = np.empty(nbins, dtype=object)
    for i in range(nbins):
        parts['logmassratios'][i] = pymc.Uninformative('logmassbias_%d' % i)
        parts['c200s'][i]= pymc.Uniform('c200_%d' % i, 1.1, 19.9)

    parts['binnedclusters'] = []
    for i in range(nbins):
        curbin = {}
        ### SET UP DATA IN CORRECT ORDERING HERE
        parts['binnedclusters'].append(curbin)

        
#####################



def buildLikelihood(parts):


    

    @pymc.observed
    def likelihood(value = parts['binnedclusters'], logmassratios = parts['logmassratios'], c200s = parts['c200s']):

        




    


####################


def runChains(model):
    pass
    
    
