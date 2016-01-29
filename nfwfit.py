#!/usr/bin/env python
#######################
# Runs a bootstrapped NFW model fit to simulated shear data.
# Galaxies are binned into an average shear profile before fitting to NFW model.
# Options to explore radial fit range, mass-concentration relation, and binning scheme in fits.
########################

import cPickle, sys, os
import numpy as np
import astropy.io.fits as pyfits
import nfwutils, bashreader, ldac
import nfwmodeltools as tools
import varcontainer
import fitmodel
import pymc
import pymc_mymcmc_adapter as pma
import scipy.integrate
import profilemaker
import simutils


#######################



def buildModel(config):


    try:
        massconRelation = simutils.buildObject(config.massconmodule, config.massconrelation, config)
    except AttributeError:
        massconRelation = None

    if massconRelation is None:
        model = NFW_Model(config = config)
    else:
        model = NFW_MC_Model(massconRelation, config = config)

    return model

    

def buildFitter(config):


    model = buildModel(config)
    

    fitter = NFWFitter(model = model, config = config)

    return fitter

####

def buildSimReader(config):

   return buildObject(config.readermodule, config.readerclass, config = config)

    
        
########################

class NFW_Model(object):

    def __init__(self, config = None):

        self.massScale = 1e14
        self.overdensity = 200
        self.config = config

        if config is not None and 'massprior' in config and config.massprior != 'linear':
            self.m200_low = 1e10
            self.m200_high = 1e17
        else:
            self.m200_low = -1e16
            self.m200_high = 1e16


        self.c200_low = 0.1
        self.c200_high = 30.

    def paramLimits(self):

        return {'m200' : (self.m200_low/self.massScale,self.m200_high/self.massScale),
                'c200' : (self.c200_low, self.c200_high)}

    def guess(self):

        guess = [10**(np.random.uniform(14, 15.5)),
                 np.random.uniform(1., 20.)]

        guess[0] = guess[0] / self.massScale

        return guess


    def setData(self, beta_s, beta_s2, zcluster):
        
        self.beta_s = beta_s
        self.beta_s2 = beta_s2
        self.zcluster = zcluster
        self.rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2




    def makeMCMCModel(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster, delta = 200):

        self.setData(beta_s, beta_s2, zcluster)

        parts = {}

        if 'massprior' in self.config and self.config.massprior == 'linear':
            parts['scaledmdelta'] = pymc.Uniform('scaledmdelta', self.m200_low/self.massScale, self.m200_high/self.massScale)
            
            @pymc.deterministic(trace=True)
            def mdelta(scaledmdelta = parts['scaledmdelta']):
                return self.massScale*scaledmdelta
            parts['mdelta'] = mdelta

        else:
            parts['logMdelta'] = pymc.Uniform('logMdelta', np.log(self.m200_low), np.log(self.m200_high))
            
            @pymc.deterministic(trace=True)
            def mdelta(logMdelta = parts['logMdelta']):

                return np.exp(logMdelta)
            parts['mdelta'] = mdelta

        parts['cdelta'] = pymc.Uniform('cdelta', self.c200_low, self.c200_high)

        rho_c = nfwutils.global_cosmology.rho_crit(zcluster)
        rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster) * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2

        @pymc.observed
        def data(value = 0.,
                 r_mpc = r_mpc,
                 ghat = ghat,
                 sigma_ghat = sigma_ghat,
                 beta_s = beta_s,
                 beta_s2 = beta_s2,
                 rho_c = rho_c,
                 rho_c_over_sigma_c = rho_c_over_sigma_c,
                 mdelta = parts['mdelta'],
                 cdelta = parts['cdelta']):

            try:
            
                logp= tools.shearprofile_like(mdelta,
                                              cdelta,
                                              r_mpc,
                                              ghat,
                                              sigma_ghat,
                                              beta_s,
                                              beta_s2,
                                              rho_c,
                                              rho_c_over_sigma_c,
                                              delta)

                return logp

            except (ValueError, ZeroDivisionError):
                
                raise pymc.ZeroProbability



        parts['data'] = data

        return pymc.Model(parts)
            


    def __call__(self, x, m200, c200):

        if m200 == 0.:
            return np.zeros_like(x)

        isNegative=m200 < 0
        if isNegative:
            m200 = np.abs(m200)

        
            
        r_scale = nfwutils.rscaleConstM(m200*self.massScale, c200, self.zcluster, self.overdensity)
    
        
        nfw_shear_inf = tools.NFWShear(x, c200, r_scale, self.rho_c_over_sigma_c)
        nfw_kappa_inf = tools.NFWKappa(x, c200, r_scale, self.rho_c_over_sigma_c)

        if isNegative:
            nfw_shear_inf = -nfw_shear_inf
        
        g = self.beta_s*nfw_shear_inf / (1 - ((self.beta_s2/self.beta_s)*nfw_kappa_inf) )

        return g


###########################################################

class NFW_MC_Model(NFW_Model):

    def __init__(self, massconRelation, config = None):

        super(NFW_MC_Model, self).__init__(config = config)
        self.massconRelation = massconRelation

    def guess(self):

        guess = [10**(np.random.uniform(14, 15.5))]

        guess[0] = guess[0] / self.massScale

        return guess


    def paramLimits(self):

        limits = super(NFW_MC_Model, self).paramLimits()

        return {'m200' : limits['m200']}

    
    def makeMCMCModel(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster):

        self.setData(beta_s, beta_s2, zcluster)

        parts = {}
        
        if 'massprior' in self.config and self.config.massprior == 'linear':
            parts['scaledm200'] = pymc.Uniform('scaledm200', self.m200_low/self.massScale, self.m200_high/self.massScale, value = self.guess()[0])
            
            @pymc.deterministic(trace=True)
            def m200(scaledm200 = parts['scaledm200']):
                return self.massScale*scaledm200
            parts['m200'] = m200

        else:
            parts['logM200'] = pymc.Uniform('logM200', np.log(self.m200_low), np.log(self.m200_high))
            
            @pymc.deterministic(trace=True)
            def m200(logM200 = parts['logM200']):

                return np.exp(logM200)
            parts['m200'] = m200

        @pymc.deterministic
        def c200(m200 = parts['m200'], zcluster = zcluster):
            
            return self.massconRelation(m200*nfwutils.global_cosmology.h, self.zcluster, self.overdensity)        
        parts['c200'] = c200

        @pymc.potential
        def c200pot(c200 = parts['c200']):
            if not np.isfinite(c200):
                raise pymc.ZeroProbability
            return 0.
        parts['c200pot'] = c200pot

        rho_c = np.float64(nfwutils.global_cosmology.rho_crit(zcluster))
        rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2

        @pymc.observed
        def data(value = 0.,
                 r_mpc = r_mpc,
                 ghat = ghat,
                 sigma_ghat = sigma_ghat,
                 beta_s = beta_s,
                 beta_s2 = beta_s2,
                 rho_c = rho_c,
                 rho_c_over_sigma_c = rho_c_over_sigma_c,
                 m200 = m200,
                 c200 = c200):


            beta_s = beta_s.astype(np.float64)
            beta_s2 = beta_s2.astype(np.float64)


            
            logprob =  tools.shearprofile_like(m200,
                                          c200,
                                          r_mpc,
                                          ghat,
                                          sigma_ghat,
                                          beta_s,
                                          beta_s2,
                                          rho_c,
                                          rho_c_over_sigma_c)


            if not np.isfinite(logprob):
                raise pymc.ZeroProbability
            return logprob

        parts['data'] = data

        return pymc.Model(parts)



    def __call__(self, x, m200):

        c200 = self.massconRelation(np.abs(m200)*self.massScale*nfwutils.global_cosmology.h, self.zcluster, self.overdensity)        

        return super(NFW_MC_Model, self).__call__(x, m200, c200)

###############################


class NFWFitter(object):

    def __init__(self, model, config):

        self.model = model
        self.config = config


    ######


    def explorePosterior(self, catalog, deltas = [200, 500, 2500]):

        r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = self.prepData(catalog)

        chains = {}

        for delta in deltas:

            mcmc_model = None
            for i in range(20):
                try:
                    mcmc_model = self.model.makeMCMCModel(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens, delta = delta)
                    break
                except pymc.ZeroProbability:
                    pass
            if mcmc_model is None:
                raise pymc.ZeroProbability

            manager = varcontainer.VarContainer()
            options = varcontainer.VarContainer()
            manager.options = options

            options.singlecore = True
            options.adapt_every = 100
            options.adapt_after = 100
            options.nsamples = 30000
            if 'nsamples' in self.config:
                options.nsamples = self.config.nsamples
            manager.model = mcmc_model

            runner = pma.MyMCMemRunner()
            runner.run(manager)
            runner.finalize(manager)

            reducedchain = dict(cdelta = np.hstack(manager.chain['cdelta'][5000::2]).astype(np.float32),
                                mdelta = np.hstack(manager.chain['mdelta'][5000::2]).astype(np.float32),
                                likelihood = np.hstack(manager.chain['likelihood'][5000::2]).astype(np.float32))

            chains[delta] = reducedchain


        return chains

        

    ######

    def minChisqMethod(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster, 
                   guess = [],
                   useSimplex=False):

        if guess == []:
            guess = self.model.guess()

        print 'GUESS: %f' % guess[0]

        self.model.setData(beta_s, beta_s2, zcluster)

        fitter = fitmodel.FitModel(r_mpc, ghat, sigma_ghat, self.model,
                                   guess = guess)
        fitter.m.limits = self.model.paramLimits()
        fitter.fit(useSimplex = useSimplex)
        if fitter.have_fit:

            uncert = fitter.uncert()
            
            return fitter.par_vals, fitter.par_err
        return None


    #######

    class BadPDFException(Exception): pass

    def scanPDF(self, catalog, config, masses = np.arange(-1.005e15, 6e15, 1e13), deltas = [200, 500, 2500]):

        if 'scanpdf_minmass' in config:
            masses = np.arange(config.scanpdf_minmass, config.scanpdf_maxmass, config.scanpdf_massstep)

        #only want to define a scan for a 1d model at this point.
        assert(isinstance(self.model, NFW_MC_Model))

        r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = self.prepData(catalog)

        self.model.setData(beta_s, beta_s2, zlens)

        fitter = fitmodel.FitModel(r_mpc, ghat, sigma_ghat, self.model,
                                   guess = self.model.guess())

        pdfs = {}

        for delta in deltas:

            chisqs = np.zeros(len(masses))

            if delta == 200:
                workingmasses = masses
            if delta != 200:
                workingmasses = np.zeros_like(masses)
                for i, curm in enumerate(masses):
                    c200 = self.model.massconRelation(np.abs(curm)*nfwutils.global_cosmology.h, 
                                                      zlens, float(delta))
                    rscale = nfwutils.rscaleConstM(np.abs(curm), c200, zlens, float(delta))
                    m200 = nfwutils.Mdelta(rscale, c200, zlens, 200)
                    if curm < 0:
                        m200 = -m200
                    workingmasses[i] = m200

            for i, mass in enumerate(workingmasses):

                chisqs[i] = fitter.statfunc(fitter.ydata,
                                                 fitter.yerr,
                                                 fitter.model(fitter.xdata,
                                                                   mass / self.model.massScale))


            pdf = np.exp(-0.5*(chisqs - np.min(chisqs)))
            pdf = pdf/scipy.integrate.trapz(pdf, masses)
            pdfs[delta] = pdf

            if np.any(np.logical_not(np.isfinite(pdf))):
                raise BadPDFException

        return (masses, pdfs)

    #######

    def runUntilNotFail(self, catalog, config):

        r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = self.prepData(catalog)

        for i in range(config.nbootstraps):

            try:
                fitresult = self.minChisqMethod(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens)


                if fitresult is not None:
                    return fitresult

            except ValueError:
                pass

        #one last try with the SIMPLEX algorithm
        return self.minChisqMethod(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens, useSimplex=True)


    #####
            

    def bootstrapFit(self, catalog, config):

        fitresults = []
        nfail = 0


        for i in range(config.nbootstraps):

            if i == 0:
                curCatalog = catalog
            else:
                curBootstrap = np.random.randint(0, len(catalog), size=len(catalog))
                curCatalog = catalog.filter(curBootstrap)

            r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = self.prepData(curCatalog, config)


            try:
                fitresult = self.minChisqMethod(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens)

                

                if fitresult is None:
                    nfail += 1
                else:
                    fitresults.append(fitresult)

            except ValueError:
                nfail += 1
                



        return fitresults, nfail




########################

def savefit(bootstrap_vals, outputname):

    with open(outputname, 'wb') as output:

        print '!!!!!!!!!!!!!!!!!!!!!!!', outputname

        cPickle.dump(bootstrap_vals, output, -1)


########################

def runNFWFit(catalogname, configname, outputname):

    try:

        config = readConfiguration(configname)

        simreader = buildSimReader(config)

        nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

        catalog = readSimCatalog(catalogname, simreader, config)

        fitter = buildFitter(config)

        if 'fitter' in config and config.fitter == 'maxlike':
            fitvals = fitter.runUntilNotFail(catalog, config)
        elif 'fitter' in config and config.fitter == 'pdf':
            fitvals = fitter.scanPDF(catalog, config)
        else:
            fitvals = fitter.explorePosterior(catalog)
    
        savefit(fitvals, outputname)

    except TypeError:
        raise TypeError(configname)

############################

class FailedCreationException(Exception): pass


if __name__ == '__main__':


    catname = sys.argv[1]
    configname = sys.argv[2]
    outname = sys.argv[3]
    

    runNFWFit(catname, configname, outname)

    if not os.path.exists(outname):
        raise FailedCreationException


