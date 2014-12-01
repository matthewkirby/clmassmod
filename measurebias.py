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

__NPROCS__ = 8
__singlecore__ = False
__samples__ = 25000

__logmass_scale__ = np.log(1e14)

#pool = Pool(__NPROCS__)


####################

def loadClusterData(answerfile, chaindir):
    # loads M-C Chains for individual clusters

    with open(answerfile, 'rb') as input:
        answers = cPickle.load(input)

    clusters = []
    for chainfile in glob.glob('%s/*.out' % chaindir):
        cluster = {}
        base = os.path.basename(chainfile)
        root, ext = os.path.splitext(base)
        
        log_mtrue = np.log(answers[root]['m200']) - __logmass_scale__

        with open(chainfile, 'rb') as chaindat:
            chain = cPickle.load(chaindat)

        cluster['id'] = root
        cluster['log_mtrue'] = log_mtrue
        cluster['like_samples'] = np.column_stack([chain['logM200'] - __logmass_scale__,
                                                   chain['c200']])

        #priors used in chain sample had flat linear c200 prior, log m200 priors
        #and our model is in logc200
        cluster['weights'] = 1./chain['c200']
        cluster['weights'] = cluster['weights'] / np.sum(cluster['weights'])


        clusters.append(cluster)

    return clusters


    

####################


def createMassBinModel(clusters, parts = None, massbinedges = np.logspace(np.log10(1e14), np.log10(5e15), 10)):
    #constant spaced log mass bins

    if parts is None:
        parts = {}

    parts['clusters'] = clusters

    nbins = len(massbinedges) - 1

    parts['bin_logmassratios'] = np.empty(nbins, dtype=object)
    parts['bin_c200s'] = np.empty(nbins, dtype=object)
    parts['bin_mass_scatter'] = np.empty(nbins, dtype=object)
    parts['bin_c200_scatter'] = np.empty(nbins, dtype=object)
    parts['bin_mc_covar'] = np.empty(nbins, dtype=object)
    for i in range(nbins):
        parts['bin_logmassratios'][i] = pymc.Uninformative('bin_logmassratio_%d' % i, value = np.log(np.random.uniform(0.5, 1.5)))
        parts['bin_c200s'][i]= pymc.Uniform('bin_c200_%d' % i, 1.1, 19.9, value = np.random.uniform(2., 6.))
        parts['bin_mass_scatter'][i] = pymc.Uniform('bin_mass_scatter_%d' % i, 0.05, 1., value=np.random.uniform(0.2, 0.5))
        parts['bin_c200_scatter'][i] = pymc.Uniform('bin_c200_scatter_%d' % i, 0.05, 1., value=np.random.uniform(0.2, 0.5))
        parts['bin_mc_covar'][i] = pymc.Uniform('bin_mc_covar_%d' % i, -1., 1., value=np.random.uniform(0.2, 0.5))


    parts['bin_assignment'] = -1*np.ones(len(clusters), dtype=np.int_)
    m_trues = np.exp(np.array([cluster['log_mtrue'] for cluster in clusters]))
    for i in range(nbins):
        selection = np.logical_and(massbinedges[i] <= m_trues, m_trues < massbinedges[i])
        parts['bin_assignment'][selection] = i

    parts['clusters'] = [clusters[i] for i in range(len(clusters)) if parts['bin_assignment'][i] != -1]
    parts['bin_assignment'] = parts['bin_assignment'][parts['bin_assignment'] != -1]

    @pymc.observed
    def clusterlikelihood(value = 0.,
                   clusters = parts['clusters'],
                   bin_assignment = parts['bin_assignment'],
                   bin_logmassratios = parts['bin_logmassratios'],
                   bin_c200s = parts['bin_c200s'],
                   bin_mass_scatter = parts['bin_mass_scatter'],
                   bin_c200_scatter = parts['bin_c200_scatter'],
                   bin_mc_covar = parts['bin_mc_covar']):

        nbins = len(bin_c200s)
        nclusters = len(clusters)

        bin_covars = [np.array([[bin_mass_scatter[i]**2, bin_mc_covar[i]*bin_mass_scatter[i]*bin_c200_scatter[i]],
                                [bin_mc_covar[i]*bin_mass_scatter[i]*bin_c200_scatter[i], bin_c200_scatter[i]**2]]) \
                          for i in range(nbins)]

        bin_invcovars = [np.linalg.inv(bin_covars[i]) for i in range(nbins)]
        bin_sqrtdetcovars = [np.sqrt(np.linalg.det(bin_covars[i])) for i in range(nbins)]

        arglist= [dict(cluster_like_samples = clusters[i]['like_samples'],
                       cluster_mean = np.array([bin_logmassratios[bin_assignment[i]] + clusters[i]['log_mtrue'],
                                                np.log(bin_c200s[bin_assignment[i]])]),
                       cluster_invcovar = bin_invcovars[bin_assignment[i]],
                       cluster_detcovar = bin_sqrtdetcovars[bin_assignment[i]]) \
                      for i in range(nclusters)]

#        cluster_logprobs = np.array(pool.map(measurebiashelper.LogSum2DGaussianWrapper,arglist))
        cluster_logprobs = np.array(map(measurebiashelper.LogSum2DGaussianWrapper,arglist))

        return np.sum(cluster_logprobs)
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
