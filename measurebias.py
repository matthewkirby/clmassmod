#!/usr/bin/env python
####################
# Python MCMC program that uses pre-sampled M-C probabilities for individual clusters 
#  to measure the bias and M-C relation for lensing
#
#  Main approach is to importance sample the integrals needed, and to use the multiprocessing module for speed
####################

import cPickle, glob, os, re, sys
import numpy as np
import pymc
import measurebiashelper
import varcontainer
import pymc_mymcmc_adapter as pma
from multiprocessing import Pool

#####################

__NPROCS__ = 4
__singlecore__ = False
__samples__ = 30000

__logmass_scale__ = np.log(1e14)


pool = None




####################

def loadClusterData(answerfile, chaindir):
    # loads M-C Chains for individual clusters

    with open(answerfile, 'rb') as input:
        answers = cPickle.load(input)

    clusters = []
    for chainfile in glob.glob('%s/*.out' % chaindir):
        print 'Loading ', chainfile
        cluster = {}
        base = os.path.basename(chainfile)
        root, ext = os.path.splitext(base)
        root, ext = os.path.splitext(root)
        
        log_mtrue = np.log(answers[root]['m200']) - __logmass_scale__

        with open(chainfile, 'rb') as chaindat:
            chain = cPickle.load(chaindat)

        logM200samples = np.hstack(chain['logM200'])[5000::10]
        logC200samples = np.log(np.hstack(chain['c200'])[5000::10])

        cluster['id'] = root
        cluster['log_mtrue'] = log_mtrue
        cluster['logmass_samples'] = logM200samples - __logmass_scale__
        cluster['logc200_samples'] = logC200samples

        #priors used in chain sample had flat linear c200 prior, log m200 priors
        #and our model is in logc200
        cluster['weights'] = 1./logC200samples
        cluster['weights'] = cluster['weights'] / np.sum(cluster['weights'])


        clusters.append(cluster)

    return clusters


    

####################




def createMassBinModel(clusters, parts = None, massbinedges = np.logspace(np.log10(1e14), np.log10(5e15), 10)):
    #constant spaced log mass bins

    massbinedges = np.log(massbinedges) - __logmass_scale__ #bring mass bins onto normalized scale

    if parts is None:
        parts = {}

    parts['clusters'] = clusters

    nbins = len(massbinedges) - 1

    parts['bin_logmassratios'] = np.empty(nbins, dtype=object)
    parts['bin_c200s'] = np.empty(nbins, dtype=object)
    parts['bin_mass_logscatter'] = np.empty(nbins, dtype=object)
    parts['bin_c200_logscatter'] = np.empty(nbins, dtype=object)
    parts['bin_mc_covar'] = np.empty(nbins, dtype=object)
    for i in range(nbins):
        parts['bin_logmassratios'][i] = pymc.Uninformative('bin_logmassratio_%d' % i, value = np.log(np.random.uniform(0.5, 1.5)))
        parts['bin_c200s'][i]= pymc.Uniform('bin_c200_%d' % i, 1.1, 19.9, value = np.random.uniform(2., 6.))
        parts['bin_mass_logscatter'][i] = pymc.Uniform('bin_mass_logscatter_%d' % i, np.log(1e-2), np.log(1.))
        parts['bin_c200_logscatter'][i] = pymc.Uniform('bin_c200_logscatter_%d' % i, np.log(1e-2), np.log(1.))
        parts['bin_mc_covar'][i] = pymc.Uniform('bin_mc_covar_%d' % i, -1., 1.)


    bin_assignment = -1*np.ones(len(clusters), dtype=np.int_)
    log_m_trues = np.array([cluster['log_mtrue'] for cluster in clusters])
    for i in range(nbins):
        selection = np.logical_and(massbinedges[i] <= log_m_trues, log_m_trues < massbinedges[i+1])
        bin_assignment[selection] = i
    parts['bin_assignment'] = bin_assignment
    
    

    parts['clusters'] = [clusters[i] for i in range(len(clusters)) if parts['bin_assignment'][i] != -1]
    parts['bin_assignment'] = parts['bin_assignment'][parts['bin_assignment'] != -1]
    parts['clusters_inbin'] = [np.arange(len(parts['clusters']))[parts['bin_assignment'] == i] for i in range(nbins)]
    

    measurebiashelper.datastore.log_mtrues = [cluster['log_mtrue'] for cluster in parts['clusters']]
    measurebiashelper.datastore.logmass_samples = [cluster['logmass_samples'] for cluster in parts['clusters']]
    measurebiashelper.datastore.logc200_samples = [cluster['logc200_samples'] for cluster in parts['clusters']]
    measurebiashelper.datastore.weights = [cluster['weights'] for cluster in parts['clusters']]
    measurebiashelper.datastore.bin_assignments = parts['bin_assignment']
    measurebiashelper.datastore.clusters_inbin = parts['clusters_inbin']

    ##pool support
    pool = Pool(__NPROCS__)


    @pymc.observed
    def clusterlikelihood(value = 0.,
                   bin_logmassratios = parts['bin_logmassratios'],
                   bin_c200s = parts['bin_c200s'],
                   bin_mass_logscatter = parts['bin_mass_logscatter'],
                   bin_c200_logscatter = parts['bin_c200_logscatter'],
                   bin_mc_covar = parts['bin_mc_covar']):

        nbins = len(bin_c200s)

        args = [dict(number = i,
                     logmassratio = bin_logmassratios[i],
                     c200 = bin_c200s[i],
                     mass_scatter = np.exp(bin_mass_logscatter[i]),
                     c200_scatter = np.exp(bin_c200_logscatter[i]),
                     mc_covar = bin_mc_covar[i]) for i in range(nbins)]

#        cluster_logprob_partialsums = np.array(map(measurebiashelper.LogSum2DGaussianWrapper, args))
        cluster_logprob_partialsums = np.array(pool.map(measurebiashelper.LogSum2DGaussianWrapper, args))




        return np.sum(cluster_logprob_partialsums)
    parts['clusterlikelihood'] = clusterlikelihood

    return pymc.Model(parts)

    
        
#####################




def runSampler(model, outfile):

    manager = varcontainer.VarContainer()
    options = varcontainer.VarContainer()
    manager.options = options

    options.singlecore = __singlecore__
    options.adapt_every = 100
    options.adapt_after = 200
    options.outputFile = outfile
    options.nsamples = __samples__
    manager.model = model

    runner = pma.MyMCRunner()
    runner.run(manager)
    runner.finalize(manager)


    
    
###########################

def main(answerfile, chaindir, outfile):

    clusters = loadClusterData(answerfile, chaindir)
    
    model = None
    for i in range(10):
        try:
            model = createMassBinModel(clusters)
            break
        except pymc.ZeroProbability:
            continue

    runSampler(model, outfile)

##########################

if __name__ == '__main__':

    answerfile, chaindir, outfile = sys.argv[1:]
    main(answerfile, chaindir, outfile)
    

    if pool is not None:
        pool.close()
        pool.join()
