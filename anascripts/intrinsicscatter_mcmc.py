#######################
# Run a mcmc to constrain the ratio
#  and intrinsic scatter between two measures
####################################

import numpy as np
import pymc, pymc_mymcmc_adapter as pma
import varcontainer

import stats


#############################

def definePivot(x):

    pivot = np.mean(x)

    return pivot, x - pivot

#####

def recoverOffset(pivot, slopes, offsets):

    return offsets - slopes*pivot
    


###################################

def constmodel(parts, ndat, offset = (0., 2.0)):

    offset = pymc.Uniform('offset', offset[0], offset[1])
    parts['offset'] = offset


    @pymc.deterministic(trace = False)
    def model(b = offset):
        return b*np.ones(ndat)
    parts['model'] = model
    

################               


#parts is a dictionary
#pivot data before setting up model

def linearmodel(x, parts, thetarange = (-np.pi/2., np.pi/2.), offset = (-1., 1.), forcePositive = True):

#    theta = pymc.Uniform('theta', thetarange[0], thetarange[1])
#    parts['theta'] = theta
#    
#    @pymc.deterministic
#    def slope(theta = theta):
#        return np.tan(theta)
#    parts['slope'] = slope
#
#
    slope = pymc.Uniform('slope', -10., 10.)
    parts['slope'] = slope

    try:
        
        offset = pymc.Uniform('offset', offset[0], offset[1])

    except TypeError:

        offset = offset
    parts['offset'] = offset


    @pymc.deterministic(trace = False)
    def model(x=x, m= slope, b = offset):
        return m*x +b
    parts['model'] = model

    if forcePositive is True:
        @pymc.potential
        def positive(model = model):
            if (model <= 0).any():
                raise pymc.ZeroProbability
            return 0.
        parts['positive'] = positive

    
###########################        




def intrinsicScatter_gausserr(x, y, yerr, parts, sigmarange = (0.005, 1.0)):

    #x and y are dictionaries with common keys.
    # the dictionary values are numpy arrays of the values for x[key] and y[key]
    # (for example, x[key] is the lensing mass and y[key] is the SZ mass for cluser "key"

              

    sigma = pymc.Uniform('sigma', sigmarange[0], sigmarange[1])
    parts['sigma'] = sigma


    @pymc.stochastic(observed=True)
    def data(value = 1., x = x, y = y, yerr = yerr, model = parts['model'], sigma = sigma):

        convolvedtau = 1./(yerr**2 + sigma**2)

        return pymc.normal_like(x = y, mu = model, tau = convolvedtau)


    parts['data'] = data


########################

def lognormIS_chainerr(x, y, parts, logsigmarange = (np.log(1e-4), np.log(1.0))):

    logsigma = pymc.Uniform('logsigma', logsigmarange[0], logsigmarange[1])
    parts['logsigma'] = logsigma

    logy = [np.log(yi) for yi in y]

    @pymc.deterministic(trace = False)
    def logmodel(model = parts['model']):
        try:
            return np.log(model)
        except FloatingPointError, e:
            print model
            raise e
    parts['logmodel'] = logmodel

    print [np.mean(yi) for yi in y]

    @pymc.stochastic(observed=True)
    def data(value = 1., x = x, logy = logy, logmodel = logmodel, logsigma = logsigma):

        sigma = np.exp(logsigma)

        logp = 0.
        for i in range(len(x)):
            logp += stats.LogSumGaussian(logy[i], logmodel[i], sigma)

        return logp


    parts['data'] = data

#########################

def const_lognormIS_gaussx_chainy(y, parts, logsigmarange = (np.log(1e-2), np.log(1.0))):
    #call const first

    logsigma = pymc.Uniform('logsigma', logsigmarange[0], logsigmarange[1])
    parts['logsigma'] = logsigma

    logy = [np.log(yi) for yi in y]

    @pymc.deterministic(trace = False)
    def logmodel(model = parts['model']):
        try:
            return np.log(model)
        except FloatingPointError, e:
            print model
            raise e
    parts['logmodel'] = logmodel


    @pymc.stochastic(observed=True)
    def data(value = 1., logy = logy, logmodel = logmodel, logsigma = logsigma):

        sigma = np.exp(logsigma)

        logp = 0.
        for i in range(len(y)):
            logp += stats.LogSumGaussian(logy[i], logmodel, sigma)

        return logp


    parts['data'] = data




#########################

def linear_lognormIS_gaussx_chainy(x,xerr,y, parts, logsigmarange = (np.log(1e-2), np.log(1.0))):
    #call linear first

    logsigma = pymc.Uniform('logsigma', logsigmarange[0], logsigmarange[1])
    parts['logsigma'] = logsigma

    logy = [np.log(yi) for yi in y]

    @pymc.deterministic(trace = False)
    def logmodel(model = parts['model']):
        try:
            return np.log(model)
        except FloatingPointError, e:
            print model
            raise e
    parts['logmodel'] = logmodel


    @pymc.stochastic(observed=True)
    def data(value = 1., x = x, xerr = xerr, logy = logy, logmodel = logmodel, slope = parts['slope'], logsigma = logsigma):

        sigma_int = np.exp(logsigma)
        sigma = np.sqrt(sigma_int**2 + (slope*xerr/x)**2)

        logp = 0.
        for i in range(len(x)):
            logp += stats.LogSumGaussian(logy[i], logmodel[i], sigma[i])

        return logp


    parts['data'] = data

#########################

def linear_normIS_gaussx_chainy(x,xerr,y, parts, logsigmarange = (np.log(1e-2), np.log(1.0))):
    #call linear first

    logsigma = pymc.Uniform('logsigma', logsigmarange[0], logsigmarange[1])
    parts['logsigma'] = logsigma


    @pymc.stochastic(observed=True)
    def data(value = 1., x = x, xerr = xerr, y = y, model = parts['model'], slope = parts['slope'], logsigma = logsigma):

        sigma_int = np.exp(logsigma)
        sigma = np.sqrt(sigma_int**2 + (slope*xerr)**2)

        logp = 0.
        for i in range(len(x)):
            logp += stats.LogSumGaussian(y[i], model[i], sigma[i])

        return logp


    parts['data'] = data



#########################

def KellyModel(x, xerr, y, yerr, xycovar, parts, ngauss = 3):
    #Implementation of Kelly07 model, but without nondetection support



    #Prior as defined in section 6.1 of Kelly07

    alpha = pymc.Uninformative('alpha', value = np.random.uniform(-1, 1))
    parts['alpha'] = alpha
    beta = pymc.Uninformative('beta', value = np.random.uniform(-np.pi/2, np.pi/2))
    parts['beta'] = beta

    sigint2 = pymc.Uniform('sigint2', 1e-4, 1.)
    parts['sigint2'] = sigint2
    
    piprior = pymc.Dirichlet('pi', np.ones(ngauss))
    parts['piprior'] = piprior

    @pymc.deterministic(trace=False)
    def pis(piprior = piprior):
        lastpi = 1. - np.sum(piprior)
        allpi = np.zeros(ngauss)
        allpi[:-1] = piprior
        allpi[-1] = lastpi
        return allpi
    parts['pis'] = pis

    mu0 = pymc.Uninformative('mu0', np.random.uniform(-1, 1))
    parts['mu0'] = mu0
    w2 = pymc.Uniform('w2', 1e-4, 1e4)
    parts['w2'] = w2


    xvars = pymc.InverseGamma('xvars', 0.5, w2, size=ngauss+1)  #dropping the 1/2 factor on w2, because I don't think it matters
    parts['xvars'] = xvars

    @pymc.deterministic(trace=False)
    def tauU2(xvars = xvars):
        return 1./xvars[-1]
    parts['tauU2'] = tauU2

    xmus = pymc.Normal('xmus', mu0, tauU2, size = ngauss)
    parts['xmus'] = xmus

    

    @pymc.observed
    def likelihood(value = 0., x = x, xerr2 = xerr**2, y = y, yerr2 = yerr**2, xycovar = xycovar, 
                   alpha = alpha, beta = beta, sigint2 = sigint2, pis = pis, xmus = xmus, xvars = xvars):

        return stats.kelly_like(x = x, 
                                xerr2 = xerr2,
                                y = y,
                                yerr2 = yerr2,
                                xycovar = xycovar,
                                alpha = alpha,
                                beta = beta,
                                sigint2 = sigint2,
                                pis = pis,
                                mus = xmus,
                                tau2 = xvars[:-1])

    parts['likelihood'] = likelihood


    
    
    



    

    

    
    


#########################


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
    
