#!/usr/bin/env python
#######################
# Runs an NFW model fit to simulated shear data.
# Galaxies are binned into an average shear profile before fitting to NFW model.
# Options to explore radial fit range, mass-concentration relation, and binning scheme in fits.
########################

import cPickle, sys, os
import numpy as np
import astropy.io.fits as pyfits
import nfwutils, bashreader, ldac
import nfwmodeltools as tools
import varcontainer
import pymc
import pymc_mymcmc_adapter as pma
import scipy.integrate
import profilebuilder
import simutils


#######################


class BasicNFWBiasModel(object):

    def configure(self, config):

        self.generate_nclusters = config['generate_nclusters']
        self.generate_minmass = config['generate_minmass']
        self.generate_maxmass = config['generate_maxmass']
        self.generate_mcrelation = config['generate_mcrelation']

    def generateDataset(self):

        

    def buildModel(self, data):

        nclusters = len(data)

        parts = {}


        ########
        #mass bias population priors
        parts['bias_logmu'] = pymc.Uniform('bias_logmu', -1., 1.)
        parts['bias_logsigma'] = pymc.Uniform('bias_logsigma', np.log(0.05), np.log(10))

        @pymc.deterministic(trace=False)
        def bias_tau(logsigma = parts['bias_logsigma']):
            return 1./np.exp(logsigma)**2
        parts['bias_tau'] = sigma
        

        for cluster_idx in range(nclusters):
            parts['log_m200s'][cluster_idx] = np.Normal('log_m200s_%d' % cluster_idx,
                                                           np.log(data['truth'][cluster_idx]) + parts['bias_logmu'],
                                                        parts['bias_tau'])

        @pymc.deterministic(trace=False)
        def m200s(log_m200s = parts['log_m200s']):
            return np.exp(log_m200s)
        parts['m200s'] = m200s

        
    ########
    # concentration population priors
    parts['concen_logmu'] = pymc.Uniform('concen_logmu', 0., np.log(25.))
    parts['concen_logsigma'] = pymc.Uniform('concen_logsigma', np.log(0.05), np.log(10))

    @pymc.deterministic(trace=True)
    def concen_sigma(logsigma = parts['concen_logsigma']):
        return np.exp(logsigma)
    
    parts['concen_sigma'] = sigma
    
    parts['logc200s'] = np.empty(nclusters, type=object)
    for cluster_idx in range(nclusters):
        parts['logc200s'][cluster_idx] = pymc.Normal('logc200_%d' % cluster_idx,
                                                     parts['concen_logmu'],
                                                     1./parts['concen_sigma']**2)


    #########
    # likelihood

    @pymc.observed
    def likelihood(value = 0., data = data,
                   m200s = parts['m200s'],
                   logc200s = parts['logc200s']):

        logp = np.sum(map(lambda params: nfwmodeltools.shearprofile_like(params[0],  #m200
                                                                          np.exp(params[1]),  #c200
                                                                          params[2].r_mpc,
                                                                          params[2].shear,
                                                                          params[2].shearerr,
                                                                          params[2].avebeta,
                                                                          params[2].avebeta2,
                                                                          params[2].rho_c,
                                                                          params[2].rho_c_over_sigma_c,
                                                                          200.),  #delta
                           zip(m200s, log200s, data)))

        
        return logp
    parts['likelihood'] = likelihood


    return parts
        
            


####################################

 

def savefit(bootstrap_vals, outputname):

    with open(outputname, 'wb') as output:

        print '!!!!!!!!!!!!!!!!!!!!!!!', outputname

        cPickle.dump(bootstrap_vals, output, -1)


########################

def runNFWFit(catalogname, configname, outputname):

    config, simreader = preloadNFWFit(configname)

    runNFWFit_Preloaded(simreader, catalogname, config, outputname)

##########################

def preloadNFWFit(configname):

    config = simutils.readConfiguration(configname)
    simreader = config['simreader']

    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    return config, simreader

###########################

def runNFWFit_Preloaded(simreader, catalogname, config, outputname):

    sim = simreader.load(catalogname)

    profilebuilder = config['profilebuilder']
    fitter = config['fitter']

    profile = profilebuilder(sim)

    fitvals = fitter(profile)

    savefit(fitvals, outputname)

############################


class FailedCreationException(Exception): pass


if __name__ == '__main__':


    catname = sys.argv[1]
    configname = sys.argv[2]
    outname = sys.argv[3]
    

    runNFWFit(catname, configname, outname)

    if not os.path.exists(outname):
        raise FailedCreationException


