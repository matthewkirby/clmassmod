#######################
# Fit a model which is lognorm intrinsic scatter and linear norm noise
########################

import numpy as np
import pymc
import deconvolvedlognormtools as dlntools
import varcontainer
import pymc_mymcmc_adapter as pma

########################


def buildModel(mlens, merr, mtrue):

    massnorm = 1e15

    parts = {}

    parts['logmu'] = pymc.Uniform('logmu', -1., 1.)
    parts['logsigma'] = pymc.Uniform('logsigma', np.log(1e-4), np.log(10))

    @pymc.deterministic(trace=True)
    def sigma(logsigma = parts['logsigma']):
        return np.exp(logsigma)
    
    parts['sigma'] = sigma

    @pymc.observed
    def data(value = 0., mlens = mlens/massnorm, merr = merr/massnorm, mtrue = mtrue/massnorm, logmu = parts['logmu'], sigma = parts['sigma']):

        return dlntools.loglinearlike(mlens = mlens,
                                      merr = merr,
                                      mtrue = mtrue,
                                      logmu = logmu,
                                      sigma = sigma)
    parts['data'] = data

    return parts

#######################################


def runFit(parts):

    N = pymc.NormApprox(parts)
    N.fit()
    
    return (N.mu[N.logmu], N.C[N.logmu]), (N.mu[N.logsigma], N.C[N.logsigma])


#######################################

def sample(parts, outputfile, samples, adaptevery = 100, adaptafter = 100, singlecore = False):


    options = varcontainer.VarContainer()
    options.outputFile = outputfile
    options.singlecore = singlecore
    options.nsamples = samples
    options.adapt_every = adaptevery
    options.adapt_after = adaptafter

    manager = varcontainer.VarContainer()
    manager.options = options
    manager.model = pymc.Model(parts)

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
    
