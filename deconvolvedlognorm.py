#######################
# Fit a model which is lognorm intrinsic scatter and linear norm noise
########################

import numpy as np
import pymc
import deconvolvedlognormtools as dlntools
import varcontainer
import pymc_mymcmc_adapter as pma
import scipy.stats

########################


def buildModel(mlens, merr, mtrue, nsamples = 150, noisefloor = False):

    massnorm = 1e15

    if noisefloor is True:
        merr = np.maximum(np.sqrt(merr**2 + (0.2*mtrue)**2), 0.1*np.abs(mlens))

    parts = {}

    parts['logmu'] = pymc.Uniform('logmu', -1., 1.)
    parts['logsigma'] = pymc.Uniform('logsigma', np.log(1e-4), np.log(10))

    @pymc.deterministic(trace=True)
    def sigma(logsigma = parts['logsigma']):
        return np.exp(logsigma)
    
    parts['sigma'] = sigma

    nclusters = len(mlens)
    ml_ints = np.zeros((nclusters, nsamples))
    delta_logmls = np.zeros((nclusters, nsamples))
    for i in range(nclusters):
        a, b = (0. - mlens[i]) /merr[i], (1e17 - mlens[i]) / merr[i]
        ml_ints[i,:] = mlens[i] + merr[i]*scipy.stats.truncnorm.rvs(a,b,size=nsamples)
        delta_logmls[i,:] = np.log(ml_ints[i,:]) - np.log(mtrue[i])

    

    @pymc.observed
    def data(value = 0., ml_ints = ml_ints/massnorm, delta_logmls = delta_logmls, 
             logmu = parts['logmu'], sigma = parts['sigma']):

        return dlntools.loglinearlike(ml_ints = ml_ints,
                                      delta_logmls = delta_logmls,
                                      logmu = logmu,
                                      sigma = sigma)
    parts['data'] = data

    return parts


#######################################

def buildOutlierModel(mlens, merr, mtrue, nsamples = 150, outlierfactor=100):

    massnorm = 1e15

    parts = {}

    parts['logmu'] = pymc.Uniform('logmu', -1., 1.)
    parts['logsigma'] = pymc.Uniform('logsigma', np.log(1e-4), np.log(10))

    @pymc.deterministic(trace=True)
    def sigma(logsigma = parts['logsigma']):
        return np.exp(logsigma)
    
    parts['sigma'] = sigma

    parts['fracoutliers'] = pymc.Uniform('fracoutliers', 0., 1.)

    nclusters = len(mlens)
    ml_ints = np.zeros((nclusters, nsamples))
    delta_logmls = np.zeros((nclusters, nsamples))
    outlier_ml_ints = np.zeros((nclusters, nsamples))
    outlier_delta_logmls = np.zeros((nclusters, nsamples))
    for i in range(nclusters):
        a, b = (0. - mlens[i]) /merr[i], (1e17 - mlens[i]) / merr[i]
        ml_ints[i,:] = mlens[i] + merr[i]*scipy.stats.truncnorm.rvs(a,b,size=nsamples)
        delta_logmls[i,:] = np.log(ml_ints[i,:]) - np.log(mtrue[i])

        a, b = (0. - mlens[i]) /outlierfactor*merr[i], (1e17 - mlens[i]) / outlierfactor*merr[i]
        outlier_ml_ints[i,:] = mlens[i] + outlierfactor*merr[i]*scipy.stats.truncnorm.rvs(a,b,size=nsamples)
        outlier_delta_logmls[i,:] = np.log(outlier_ml_ints[i,:]) - np.log(mtrue[i])

    

    @pymc.observed
    def data(value = 0., ml_ints = ml_ints/massnorm, delta_logmls = delta_logmls,
             outlier_ml_ints = outlier_ml_ints, outlier_delta_logmls = outlier_delta_logmls,
             logmu = parts['logmu'], sigma = parts['sigma'],
             fracoutliers = parts['fracoutliers']):

        return dlntools.outlierloglinearlike(ml_ints,
                                             delta_logmls,
                                             outlier_ml_ints,
                                             outlier_delta_logmls,
                                             logmu,
                                             sigma,
                                             fracoutliers)
    parts['data'] = data

    return parts


    

    

#######################################


def runFit(parts):

    N = pymc.NormApprox(parts)
    N.fit()
    
    return (N.mu[N.logmu], np.sqrt(N.C[N.logmu])), (N.mu[N.logsigma], np.sqrt(N.C[N.logsigma]))


#######################################

def sample(parts, outputfile, samples, adaptevery = 100, adaptafter = 100, singlecore = False):


    options = varcontainer.VarContainer()
    options.outputFile = outputfile
    options.singlecore = singlecore
    options.nsamples = samples
    options.adapt_every = adaptevery
    options.adapt_after = adaptafter
    options.restore=True

    manager = varcontainer.VarContainer()
    manager.options = options
    manager.model = pymc.Model(parts)
    
    assert(np.isfinite(manager.model.logp))

    runner = pma.MyMCRunner()

    runner.run(manager)

    runner.finalize(manager)




##############################

def memsample(parts, samples, adaptevery = 100, adaptafter = 100):

    options = varcontainer.VarContainer()
    options.singlecore = True
    options.nsamples = samples
    options.adapt_every = adaptevery
    options.adapt_after = adaptafter

    manager = varcontainer.VarContainer()
    manager.options = options
    manager.model = pymc.Model(parts)

    runner = pma.MyMCMemRunner()

    runner.run(manager)

    runner.finalize(manager)

    return manager.chain


#################

