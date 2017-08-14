import pymc
import numpy as np
import concentrationfittools as cfittools


def buildMCMCModel_mc(halos, maxsamples = 2000):
    '''
    This now runs the MCMC to fit both mass and concentration

    halos is a list of dictionaries corresponding to each halo
    
    id
    true_mass
    mass_samples
    concentration_samples

    These are filled in deconvolvedlognorm methods.
    '''
    massnorm = 1e15

    parts = {}

    parts['logmu'] = pymc.Uniform('logmu', -1., 1.)  # Tracking bias
    parts['logsigma'] = pymc.Uniform('logsigma', np.log(1e-2), np.log(10))  

    parts['logmu_c'] = pymc.Uniform('logmu_c', 0., np.log(15.))  # concentrations run between 1 - 15 for clusters
    parts['logsigma_c'] = pymc.Uniform('logsigma_c', np.log(1e-2), np.log(.5)) # realistic scatter in terms of percentage
    
    @pymc.deterministic(trace=True)
    def sigma(logsigma = parts['logsigma']):
        return np.exp(logsigma)

    @pymc.deterministic(trace=True)
    def sigma_c(logsigma_c = parts['logsigma_c']) :
        return np.exp(logsigma_c)
    
    parts['sigma'] = sigma
    parts['sigma_c'] = sigma_c

    nclusters = len(halos)
    ml_ints = np.zeros((nclusters, maxsamples))
    delta_logmls = np.zeros((nclusters, maxsamples))
    ngoodsamples = np.zeros(nclusters, dtype=np.int)

    cl_ints = np.zeros((nclusters, maxsamples))
    
    for i in range(nclusters):
        mass_samples = halos[i]['mass_samples']
        concentration_samples = halos[i]['concentration_samples']
        positive_samples = mass_samples[mass_samples > 0]
        positive_samples_c = concentration_samples[mass_samples > 0] 

        navailablesamples = len(positive_samples)
        takesamples = min(navailablesamples, maxsamples)

        if navailablesamples < 25:
            print 'Need more samples: ', halos[i]['id'], navailablesamples

        # Random downsampling
        sample_indices = np.random.permutation(np.arange(navailablesamples))[:takesamples]

        ml_ints[i,:takesamples] = positive_samples[sample_indices]
        cl_ints[i,:takesamples] = positive_samples_c[sample_indices]

        # Bias
        delta_logmls[i,:takesamples] = np.log(ml_ints[i,:takesamples]) - np.log(halos[i]['true_mass'])
        ngoodsamples[i] = takesamples

    

    @pymc.observed
    def data(value = 0., ml_ints = ml_ints/massnorm, delta_logmls = delta_logmls, ngoodsamples = ngoodsamples,
             logmu = parts['logmu'], sigma = parts['sigma'],
             cl_ints = cl_ints, logmu_c = parts['logmu_c'], sigma_c = parts['sigma_c']
    ):

        return cfittools.mcmcloglinearlike_mc(ml_ints = ml_ints,
                                           delta_logmls = delta_logmls,
                                           ngoodsamples = ngoodsamples,
                                           logmu = logmu,
                                           sigma = sigma,
                                           cl_ints = cl_ints,
                                           logmu_c = logmu_c,
                                           sigma_c = sigma_c
        )
    parts['data'] = data

    return parts


