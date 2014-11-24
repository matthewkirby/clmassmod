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


def buildModel(clusters):
    #buid pymc model that infers bias, scatter, and m-c model

    parts = {}

    #start with a linear model for \vec alpha, relating log(m_true) and log(m_lens)

    parts['massnorm'] = 5e14
    parts['log_massnorm'] = np.log(massnorm)

    parts['theta'] = pymc.Uniform('theta', -np.pi/2., np.pi/2.)
    
    @pymc.deterministic
    def slope(theta = parts['theta']):
        return np.tan(theta)
    parts['slope'] = slope

    parts['offset'] = pymc.Uniform('offset', -100., 100.)

    parts['log_mtrue_norm'] = np.array([cluster['log_mtrue'] for cluster in clusters]) - parts['log_massnorm']

    @pymc.deterministic(trace = False)
    def log_mlensing_pred(x=parts['log_mtrue_norm'], m= parts['slope'], b = parts['offset']):
        return m*x +b
    parts['pred_log_mlens_norm'] = model
    
    parts['pred_c200'] = 4.*np.ones(len(clusters))

    parts['logmlens_log_scatter'] = pymc.Uniform('logmlens_log_scatter', np.log(1e-4), np.log(1.))
    @pymc.deterministic
    def logmlens_scatter(logscatter = parts['logmlens_log_scatter']):
        return np.exp(logscatter)
    parts['logmlens_scatter'] = logmlens_scatter

    parts['log_c200_scatter'] = 0.116*np.log(10)


    parts['like_log_mlens_norm'] = [cluster['log_mlens'] - parts['log_massnorm'] for cluster in clusters]
    parts['like_c200'] = [cluster['c200'] for cluster in clusters]
    #no weights are needed since I'll compute P(log_mass) and the chain used P(log_mass) = const

    

    
    


    @pymc.observed




    


####################


def runChains(model):
    pass
    
    
