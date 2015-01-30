#!/usr/bin/env python
#######################
# Runs a bootstrapped NFW model fit to simulated shear data.
# Galaxies are binned into an average shear profile before fitting to NFW model.
# Options to explore radial fit range, mass-concentration relation, and binning scheme in fits.
########################

import importlib, cPickle, sys, os
import numpy as np
import astropy.io.fits as pyfits
import astropy.io.ascii as asciireader
import nfwutils, bashreader, ldac
import nfwmodeltools as tools
import varcontainer
import fitmodel
import pymc
import pymc_mymcmc_adapter as pma


#######################

def applyMask(x_arcmin, y_arcmin, zcluster, config):

    if 'targetz' in config:
        #adjust for different redshift
        curr_angdist = nfwutils.global_cosmology.angulardist(zcluster)
        newangdist = nfwutils.global_cosmology.angulardist(config.targetz)
        ratio = curr_angdist/newangdist
        x_arcmin = x_arcmin*ratio
        y_arcmin = y_arcmin*ratio


    def squaremask(x=0,y=0,theta=0, sidelength = 1.):
        #x,y in arcminutes

        theta_rad = np.pi*theta/180.
        
        dX = x_arcmin - x   
        dY = y_arcmin - y
        
        rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
        rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

        return (np.abs(rdX) + np.abs(rdY)) <= (np.sqrt(2)*sidelength/2.)

    def rectanglemask(x=0,y=0,theta=0, xlength=1., ylength=1.):

        #x,y in arcminutes

        theta_rad = np.pi*theta/180.
        
        dX = x_arcmin - x   
        dY = y_arcmin - y
        
        rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
        rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

        return np.logical_and(np.abs(rdX) < xlength/2., np.abs(rdY) < ylength/2.)


        

    def circlemask(x=0, y=0, rad=1.):
        #arcmnutes
        
        dX = x_arcmin - x
        dY = y_arcmin - y
        
        return np.sqrt(dX**2 + dY**2) < rad

    acsmask = lambda x,y: squaremask(x,y,sidelength=3.2)
    wfc3mask = lambda x,y: squaremask(x,y,theta=45.,sidelength=(3.2*4./5.))

    acscentered = lambda : np.logical_or(acsmask(0,0), wfc3mask(0, 6.))

    offsetpointing = lambda : np.logical_or(acsmask(0.,-3.), wfc3mask(0., 3.0))

    rotatedoffset = lambda : np.logical_or(acsmask(-3.,0.), wfc3mask(3., 0.0))

    offsetmosaic = lambda : np.logical_or(offsetpointing(), rotatedoffset())

    offset3 = lambda : np.logical_or(offsetmosaic(), acscentered())

    pisco3 = lambda: np.logical_or(np.logical_or(rectanglemask(0,0,0,8,6), rectanglemask(0,7,0,6,8)), rectanglemask(0,-7,0,6,8))

    pisco4 = lambda: np.logical_or(np.logical_or(rectanglemask(-5.5, 0, 0, 8, 6), rectanglemask(5.5, 0, 0, 8, 6)), 
                                   np.logical_or(rectanglemask(0, 5.5, 0, 6, 8), rectanglemask(0, -5.5, 0,6, 8)))

    def randomoffset():
        posangle = np.random.uniform(0, 2*np.pi)
        rotmatrix = np.array([[np.cos(posangle), -np.sin(posangle)],
                              [np.sin(posangle), np.cos(posangle)]])
        acspos = np.dot(rotmatrix, np.array([0., -3]))
        wfc3pos = np.dot(rotmatrix, np.array([0., 3]))

        print acspos, wfc3pos

        rotatedoffset = np.logical_or(acsmask(acspos[0], acspos[1]), wfc3mask(wfc3pos[0], wfc3pos[1]))
        return rotatedoffset

    randomoffsetmosaic = lambda : np.logical_or(randomoffset(), offsetpointing())

    centerandoffset = lambda : np.logical_or(acscentered(), offsetpointing())

    centerandrandoffset = lambda : np.logical_or(acscentered(), randomoffset())


    acsdiag = np.sqrt(2.)*3.2/2.
    squaremosaic = lambda : np.logical_or(np.logical_or(acsmask(0., acsdiag),acsmask(0.,-acsdiag)), 
                                          np.logical_or(acsmask(-acsdiag,0.),acsmask(acsdiag,0.)))

    maskcase = {'squaremask' : lambda : squaremask(config.maskx, config.masky, config.masktheta, config.masksidelength),
                'circlemask' : lambda : circlemask(config.maskx, config.masky, config.maskrad),
                'acsmask' : lambda : acsmask(config.maskx, config.masky),
                'wfc3mask' : lambda : wfc3mask(config.maskx, config.masky),
                'acscentered' : acscentered,
                'offsetpointing' : offsetpointing,
                'rotatedoffset' : rotatedoffset,
                'offsetmosaic' : offsetmosaic,
                'offset3' : offset3,
                'pisco3' : pisco3,
                'pisco4' : pisco4,
                'randomoffset' : randomoffset,
                'randomoffsetmosaic' : randomoffsetmosaic,
                'centerandoffset' : centerandoffset,
                'centerandrandoffset' : centerandrandoffset,
                'squaremosaic' : squaremosaic}

    mask = maskcase[config.maskname]()

    return mask
            

########################

class InsufficientGalaxiesException(Exception): pass

def applyDensityMask(x_arcmin, y_arcmin, zcluster, config):
    #assumes that the input catalog is rectalinear

    targetdensity = config.nperarcmin

    max_x = np.max(x_arcmin)
    min_x = np.min(x_arcmin)
    delta_x = max_x - min_x

    max_y = np.max(y_arcmin)
    min_y = np.min(y_arcmin)
    delta_y = max_y - min_y

    area = delta_x*delta_y

    if targetdensity == -1:
        #take all
        targetnumber = len(x_arcmin)

    elif 'targetz' in config:
        #adjust for different redshift
        curr_angdist = nfwutils.global_cosmology.angulardist(zcluster)
        newangdist = nfwutils.global_cosmology.angulardist(config.targetz)
        ratio = curr_angdist/newangdist
        newarea = area*ratio**2
        targetnumber = targetdensity*newarea
    else:

        targetnumber = targetdensity*area

    availablenumber = len(x_arcmin)

    if targetnumber > availablenumber:
        raise InsufficientGalaxiesException

    accept = float(targetnumber) / availablenumber

    randomthrow = np.random.random(len(x_arcmin))

    selected = randomthrow < accept


    return randomthrow < accept



########################

def readSimCatalog(catalogname, simreader, config):

    sim = simreader.load(catalogname)


    e1 = sim.g1
    e2 = sim.g2

    centeroffsetx = 0.
    centeroffsety = 0.
    if 'coresize' in config:
        offsetcat = asciireader.read('/vol/braid1/vol1/dapple/mxxl/mxxlsims/SPT_SN_offset.dat')
        m500 = sim.m500
        if m500 == 0:
            raise ValueError
        #choose one of the 50 closest in mass at same core radius
        matchingcoresize = offsetcat[offsetcat['coresize[arcmin]'] == config.coresize]
        print 'Distro Available: %d' % len(matchingcoresize)
        deltamass = matchingcoresize['M500c'] - m500
        closestsims = np.argsort(deltamass)
        selectedsim = closestsims[np.random.uniform(0, min(50, len(deltamass)))]  
        centeroffsetx = (matchingcoresize['peak_xpix[arcmin]'] - matchingcoresize['cluster_xpix'])[selectedsim]
        centeroffsety = (matchingcoresize['peak_ypix'] - matchingcoresize['cluster_ypix'])[selectedsim]
        print 'Pointing Offset: %f %f' % (centeroffsetx, centeroffsety)
    


    delta_x = sim.x_arcmin - centeroffsetx
    delta_y = sim.y_arcmin - centeroffsety

    dL = nfwutils.global_cosmology.angulardist(sim.zcluster)
    deltax_mpc = (delta_x * dL * np.pi)/(180.*60)
    deltay_mpc = (delta_y * dL * np.pi)/(180.*60)
    r_mpc = np.sqrt(deltax_mpc**2 + deltay_mpc**2)
    cosphi = deltax_mpc / r_mpc
    sinphi = deltay_mpc / r_mpc
    
    sin2phi = 2.0*sinphi*cosphi
    cos2phi = 2.0*cosphi*cosphi-1.0

     
   
        

    E = -(e1*cos2phi+e2*sin2phi)

    b1 =  e2
    b2 = -e1
    B = -(b1*cos2phi+b2*sin2phi)

    if 'shapenoise' in config:
        E = E + config.shapenoise*np.random.standard_normal(len(E))
        B = B + config.shapenoise*np.random.standard_normal(len(B))


    visiblemask = np.ones(len(E)) == 1.
    if 'maskname' in config:
        visiblemask = applyMask(delta_x, delta_y, sim.zcluster, config)

    densitymask = np.ones(len(E)) == 1.
    if 'nperarcmin' in config:
        densitymask = applyDensityMask(delta_x, delta_y, sim.zcluster, config)

    mask = np.logical_and(visiblemask, densitymask)

    r_arcmin = sim.r_arcmin
    r_mpc = sim.r_mpc
    redshifts = sim.redshifts
    beta_s = sim.beta_s

    cols = [pyfits.Column(name = 'r_arcmin', format = 'E', array = r_arcmin),
            pyfits.Column(name = 'r_mpc', format='E', array = r_mpc),
            pyfits.Column(name = 'ghat', format='E', array = E),
            pyfits.Column(name = 'gcross', format='E', array = B),
            pyfits.Column(name = 'z', format='E', array = redshifts),
            pyfits.Column(name = 'beta_s', format = 'E', array = beta_s),
            pyfits.Column(name = 'mask', format = 'L', array = mask)]
    catalog = ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))
    catalog.hdu.header['ZLENS'] = sim.zcluster


    return catalog


########################

def readConfiguration(configname):

    config = bashreader.parseFile(configname)
    return config
    

########################

def buildObject(modulename, classname, *args, **kwds):

    if modulename.lower() == 'none' or classname.lower() == 'none':
        return None

    aModule = importlib.import_module(modulename)
    aClass = getattr(aModule, classname)
    anObject = aClass(*args, **kwds)

    return anObject
    
#####

def buildProfileBuilder(config):

    return buildObject(config.profilemodule, config.profilebuilder, config)

#####

def buildModel(config):


    try:
        massconRelation = buildObject(config.massconmodule, config.massconrelation, config)
    except AttributeError:
        massconRelation = None

    if massconRelation is None:
        model = NFW_Model(config = config)
    else:
        model = NFW_MC_Model(massconRelation, config = config)

    return model

    

def buildFitter(config):

    profileBuilder = buildProfileBuilder(config)
    model = buildModel(config)

    

    fitter = NFWFitter(profileBuilder = profileBuilder, model = model, config = config)



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
            self.m200_low = -1e18
            self.m200_high = 1e18


        self.c200_low = 1.1
        self.c200_high = 19.9

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




    def makeMCMCModel(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster):

        self.setData(beta_s, beta_s2, zcluster)

        parts = {}
        
        if 'massprior' in self.config and self.config.massprior == 'linear':
            parts['scaledm200'] = pymc.Uniform('m200', self.m200_low/self.massScale, self.m200_high/self.massScale)
            
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

        parts['c200'] = pymc.Uniform('c200', self.c200_low, self.c200_high)

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
                 m200 = parts['m200'],
                 c200 = parts['c200']):

            try:
            
                logp= tools.shearprofile_like(m200,
                                              c200,
                                              r_mpc,
                                              ghat,
                                              sigma_ghat,
                                              beta_s,
                                              beta_s2,
                                              rho_c,
                                              rho_c_over_sigma_c)

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

    def __init__(self, profileBuilder, model, config):

        self.profileBuilder = profileBuilder
        self.model = model
        self.config = config


    ######

    def prepData(self, curCatalog):
        
        beta_s = np.mean(curCatalog['beta_s'], dtype=np.float64)
        beta_s2 = np.mean(curCatalog['beta_s']**2, dtype=np.float64)
        
        r_mpc, ghat, sigma_ghat = [x.astype(np.float64) for x in self.profileBuilder(curCatalog, self.config)]

        clean = sigma_ghat > 0

        zlens = curCatalog.hdu.header['ZLENS']

        return r_mpc[clean], ghat[clean], sigma_ghat[clean], beta_s, beta_s2, zlens


    ######

    def explorePosterior(self, catalog):

        r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = self.prepData(catalog)

        mcmc_model = None
        for i in range(10):
            try:
                mcmc_model = self.model.makeMCMCModel(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens)
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
        options.nsamples = 10000
        if 'nsamples' in self.config:
            options.nsamples = self.config.nsamples
        manager.model = mcmc_model
        
        runner = pma.MyMCMemRunner()
        runner.run(manager)
        runner.finalize(manager)

        return manager.chain
        

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

        cPickle.dump(bootstrap_vals, output, -1)


########################

def runNFWFit(catalogname, configname, outputname):

    config = readConfiguration(configname)

    simreader = buildSimReader(config)

    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    catalog = readSimCatalog(catalogname, simreader, config)

    fitter = buildFitter(config)

    if 'fitter' in config and config.fitter == 'maxlike':
        fitvals = fitter.runUntilNotFail(catalog, config)
    else:
        fitvals = fitter.explorePosterior(catalog)
    
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


