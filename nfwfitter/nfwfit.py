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


class NFW_Model(object):

    def __init__(self):

        self.massScale = 1e14
        self.overdensity = 200
        self.c200_low = 0.1
        self.c200_high = 30.
        

    def configure(self, config):

        if 'massprior' in config and config['massprior'] == 'log':
            self.massprior = 'log'
            self.m200_low = 1e10
            self.m200_high = 1e17
        else:
            self.massprior = 'linear'
            self.m200_low = 1e12
            self.m200_high = 1e16



    def paramLimits(self):

        return {'m200' : (self.m200_low/self.massScale,self.m200_high/self.massScale),
                'c200' : (self.c200_low, self.c200_high)}

    def guess(self):

        guess = [10**(np.random.uniform(14, 15.5)),
                 np.random.uniform(1., 20.)]

        guess[0] = guess[0] / self.massScale

        return guess


    def setData(self, beta_s, beta_s2, zcluster, zlens = None):
        
        #handle cases where shear signal has been rescaled to a different cluster redshift
        if zlens is None:
            zlens = zcluster
        
        self.beta_s = beta_s
        self.beta_s2 = beta_s2
        self.zcluster = zcluster
        self.zlens = zlens
        self.rho_c = nfwutils.global_cosmology.rho_crit(zcluster)
        #note the mixed usage of zlens and zcluster. Lensing properties were rescaled. Mass related properties (here H2, from rho_crit) were not.
        self.rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zlens) * nfwutils.global_cosmology.beta([1e6], zlens)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2






    def makeMCMCModel(self, profile, delta = 200):

        self.setData(profile.beta_s, profile.beta_s2, profile.zcluster, zlens = profile.zlens)

        parts = {}

        if self.massprior == 'linear':
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


        @pymc.observed
        def data(value = 0.,
                 r_mpc = profile.r_mpc,
                 ghat = profile.ghat,
                 sigma_ghat = profile.sigma_ghat,
                 beta_s = profile.beta_s,
                 beta_s2 = profile.beta_s2,
                 rho_c = self.rho_c,
                 rho_c_over_sigma_c = self.rho_c_over_sigma_c,
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

    def configure(self, config):

        super(NFW_MC_Model, self).configure(config)
        self.massconRelation = config['massconRelation']

    def guess(self):

        guess = [10**(np.random.uniform(14, 15.5))]

        guess[0] = guess[0] / self.massScale

        return guess


    def paramLimits(self):

        limits = super(NFW_MC_Model, self).paramLimits()

        return {'m200' : limits['m200']}

    
    def makeMCMCModel(self, profile, delta = 200):

        self.setData(profile.beta_s, profile.beta_s2, profile.zcluster, zlens = profile.zlens)

        parts = {}
        
        if  self.massprior == 'linear':
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


        @pymc.observed
        def data(value = 0.,
                 r_mpc = profile.r_mpc,
                 ghat = profile.ghat,
                 sigma_ghat = profile.sigma_ghat,
                 beta_s = profile.beta_s,
                 beta_s2 = profile.beta_s2,
                 rho_c = self.rho_c,
                 rho_c_over_sigma_c = self.rho_c_over_sigma_c,
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


class MCMCFitter(object):

    def __init__( self ) :
        self.output_type = 'mcmc'
    
    def configure(self, config):

        self.model = config['model']
        self.deltas = [200, 500, 2500]
        self.nsamples = 10000
        if 'nsamples' in config:
            self.nsamples = config['nsamples']

    def verifyfit(self, sim, profile, fitvals, outputname, raiseException = True):
        '''Not implemented yet for MCMCs - see PDFScanner for the version in 1-d'''

        return True

        


    def __call__(self, profile):

        chains = {}

        for delta in self.deltas:

            mcmc_model = None
            for i in range(20):
                try:
                    mcmc_model = self.model.makeMCMCModel(profile, delta = delta)
                    break
                except pymc.ZeroProbability:
                    pass
            if mcmc_model is None:
                raise pymc.ZeroProbability
            # This sets up Adam Mantz's version of an MCMC sampler for production code calculations.
            # This is stored in mymcmc_adapter.py (converts to talk with other MCMC code)
            # Imported as pma above.
            
            manager = varcontainer.VarContainer()
            options = varcontainer.VarContainer()
            manager.options = options

            options.singlecore = True
            options.adapt_every = 100
            options.adapt_after = 100
            options.nsamples = self.nsamples
            manager.model = mcmc_model

            runner = pma.MyMCMemRunner()
            runner.run(manager)
            runner.finalize(manager)

            reducedchain = dict(cdelta = np.hstack(manager.chain['cdelta'][5000::2]).astype(np.float32),
                                mdelta = np.hstack(manager.chain['mdelta'][5000::2]).astype(np.float32),
                                likelihood = np.hstack(manager.chain['likelihood'][5000::2]).astype(np.float32))

            chains[delta] = reducedchain


        return chains


##########



#######
        
class BadPDFException(Exception): pass

class PDFScanner(object):

    def __init__(self) :

        self.output_type = 'pdf'

        
    def configure(self, config):

        self.model = config['model']
        self.deltas = [200, 500, 2500]

        self.masses = np.arange(-1.005e15, 6e15, 1e13)
        if 'scanpdf_minmass' in config:
            self.masses = np.arange(config['scanpdf_minmass'], config['scanpdf_maxmass'], config['scanpdf_massstep'])

    def verifyfit(self, sim, profile, fitvals, outputname, raiseException = True):
        '''Raises FailedFitException and dumps intermediates to pkl file if verify failes'''

        masses, pdfs = fitvals

        for delta in pdfs.keys():
            maxpos = np.argmax(pdfs[delta])
            if maxpos == 0 or maxpos == (len(masses)-1):
            
                if raiseException:
                    dump(sim, profile, fitvals, outputname)
                    raise FailedFitException
                return False

        return True



            

    def __call__(self, profile):


        #only want to define a scan for a 1d model at this point.
        assert(isinstance(self.model, NFW_MC_Model))

        self.model.setData(profile.beta_s, profile.beta_s2, profile.zcluster, zlens = profile.zlens)


        pdfs = {}

        masses = self.masses

        for delta in self.deltas:

            logprob = np.zeros(len(masses))

            if delta == 200:
                workingmasses = masses
                c200s = np.array([self.model.massconRelation(np.abs(curm)*nfwutils.global_cosmology.h, 
                                                      profile.zcluster, float(delta)) for curm in workingmasses])
            elif delta != 200:
                workingmasses = np.zeros_like(masses)
                c200s = np.zeros_like(masses)
                for i, curm in enumerate(masses):
                    c200 = self.model.massconRelation(np.abs(curm)*nfwutils.global_cosmology.h, 
                                                      profile.zcluster, float(delta))
                    rscale = nfwutils.rscaleConstM(np.abs(curm), c200, profile.zcluster, float(delta))
                    m200 = nfwutils.Mdelta(rscale, c200, profile.zcluster, 200)
                    if curm < 0:
                        m200 = -m200
                    workingmasses[i] = m200
                    c200s[i] = c200

            for i in range(len(workingmasses)):
                mass = workingmasses[i]
                c200 = c200s[i]

                logprob[i] = tools.shearprofile_like(mass, c200,
                                                    profile.r_mpc,
                                                    profile.ghat,
                                                    profile.sigma_ghat,
                                                    self.model.beta_s,
                                                    self.model.beta_s2,
                                                    self.model.rho_c,
                                                    self.model.rho_c_over_sigma_c,
                                                    200.)





            pdf = np.exp(logprob - np.max(logprob))
            pdf = pdf/scipy.integrate.trapz(pdf, masses)
            pdfs[delta] = pdf

            if np.any(np.logical_not(np.isfinite(pdf))):
                raise BadPDFException

        return (masses, pdfs)

    #######

def convertLikelihoodScan(model, delta, masses, pdf200, zcluster):

    #treats input pdf as a likelihood scan & rescales axis. Does not transform like a PDF!!!

    targetmasses = np.zeros_like(masses)
    for i, curm in enumerate(masses):
        c200 = model.massconRelation(np.abs(curm)*nfwutils.global_cosmology.h, 
                                     zcluster, 200.)
        r200 = nfwutils.rdeltaConstM(np.abs(curm), zcluster, 200.)
        rscale = r200/c200
        m_target = nfwutils.Mdelta(rscale, c200, zcluster, delta)
        if curm < 0:
            m_target = -m_target
        targetmasses[i] = m_target
    

    targetpdf = np.interp(masses, targetmasses, pdf200)

    targetpdf = targetpdf / np.trapz(targetpdf, masses)

    return targetpdf

    
########################

class FailedFitException(Exception): pass

def dump(sim, profile, fitvals, outputname):
    
    with open('{}.err.pkl'.format(outputname), 'wb') as output:
        cPickle.dump((sim, profile, fitvals), output, -1)


########################

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

    fitter.verifyfit(sim, profile, fitvals, outputname)

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


