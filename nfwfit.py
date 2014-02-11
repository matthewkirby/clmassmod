#!/usr/bin/env python
#######################
# Runs a bootstrapped NFW model fit to simulated shear data.
# Galaxies are binned into an average shear profile before fitting to NFW model.
# Options to explore radial fit range, mass-concentration relation, and binning scheme in fits.
########################

import importlib, cPickle, sys, os
import numpy as np
import astropy.io.fits as pyfits
import nfwutils, bashreader, ldac
import fitmodel, nfwmodeltools as tools



#######################

def applyMask(cat, config):

    def squaremask(x=0,y=0,theta=0, sidelength = 1.):
        #x,y in arcminutes

        theta_rad = np.pi*theta/180.
        
        dX = cat['X_arcmin'] - x   #column is mislabelled
        dY = cat['Y_arcmin'] - y
        
        rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
        rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

        return (np.abs(rdX) + np.abs(rdY)) <= (np.sqrt(2)*sidelength/2.)

    def circlemask(x=0, y=0, rad=1.):
        #arcmnutes
        
        dX = cat['X_arcmin'] - x    #column is mislabelled
        dY = cat['Y_arcmin'] - y
        
        return np.sqrt(dX**2 + dY**2) < rad

    acsmask = lambda x,y: squaremask(x,y,sidelength=3.2)
    wfc3mask = lambda x,y: squaremask(x,y,theta=45.,sidelength=(3.2*4./5.))

    acscentered = lambda : np.logical_or(acsmask(0,0), wfc3mask(0, 6.))

    offsetpointing = lambda : np.logical_or(acsmask(0.,-3.), wfc3mask(0., 3.0))

    rotatedoffset = lambda : np.logical_or(acsmask(-3.,0.), wfc3mask(3., 0.0))

    offsetmosaic = lambda : np.logical_or(offsetpointing(), rotatedoffset())

    offset3 = lambda : np.logical_or(offsetmosaic(), acscentered())

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
                'randomoffset' : randomoffset,
                'randomoffsetmosaic' : randomoffsetmosaic,
                'centerandoffset' : centerandoffset,
                'centerandrandoffset' : centerandrandoffset,
                'squaremosaic' : squaremosaic}

    mask = maskcase[config.maskname]()

    return cat.filter(mask)
            



    

        


########################

def readSimCatalog(catalogname, simreader, config):

    sim = simreader.load(catalogname)

    r_arcmin = np.sqrt(sim.delta_arcmin[0]**2 + sim.delta_arcmin[1]**2)
    r_mpc = np.sqrt(sim.delta_mpc[0]**2 + sim.delta_mpc[1]**2)

    cosphi = sim.delta_mpc[0] / r_mpc
    sinphi = sim.delta_mpc[1] / r_mpc
    
    sin2phi = 2.0*sinphi*cosphi
    cos2phi = 2.0*cosphi*cosphi-1.0
    

    e1 = sim.g1
    e2 = sim.g2

    E = -(e1*cos2phi+e2*sin2phi)

    b1 =  e2
    b2 = -e1
    B = -(b1*cos2phi+b2*sin2phi)

    redshifts = sim.redshifts
    beta_s = sim.beta_s

    cols = [pyfits.Column(name = 'r_arcmin', format = 'E', array = r_arcmin),
            pyfits.Column(name = 'r_mpc', format='E', array = r_mpc),
            pyfits.Column(name = 'ghat', format='E', array = E),
            pyfits.Column(name = 'gcross', format='E', array = B),
            pyfits.Column(name = 'z', format='E', array = redshifts),
            pyfits.Column(name = 'beta_s', format = 'E', array = beta_s)]
    catalog = ldac.LDACCat(pyfits.new_table(pyfits.ColDefs(cols)))
    catalog.hdu.header.update('ZLENS', sim.zcluster)


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

def buildFitter(config):

    profileBuilder = buildObject(config.profilemodule, config.profilebuilder, config)

    try:
        massconRelation = buildObject(config.massconmodule, config.massconrelation, config)
    except AttributeError:
        massconRelation = None

    fitter = NFWFitter(profileBuilder = profileBuilder, massconRelation = massconRelation)



    return fitter

####

def buildSimReader(config):

   return buildObject(config.readermodule, config.readerclass)

    
        
########################

class NFW_Model(object):

    def __init__(self):

        self.massScale = 1e14

    def paramLimits(self):

        return {'m200' : (1e13/self.massScale,1e16/self.massScale),
                'c200' : (1., 20.)}

    def guess(self):

        guess = [10**(np.random.uniform(13, 16)),
                 np.random.uniform(1., 20.)]

        guess[0] = guess[0] / self.massScale

        return guess


    def setData(self, beta_s, beta_s2, zcluster):
        
        self.beta_s = beta_s
        self.beta_s2 = beta_s2
        self.zcluster = zcluster
        self.overdensity = 200


    def __call__(self, x, m200, c200):
                
        r_scale = nfwutils.rscaleConstM(m200*self.massScale, c200, self.zcluster, self.overdensity)
        
        nfw_shear_inf = tools.NFWShear(x, c200, r_scale, self.zcluster)
        nfw_kappa_inf = tools.NFWKappa(x, c200, r_scale, self.zcluster)
        
        g = self.beta_s*nfw_shear_inf / (1 - ((self.beta_s2/self.beta_s)*nfw_kappa_inf) )
                
        return g

#######

class NFW_MC_Model(NFW_Model):

    def __init__(self, massconRelation):

        super(NFW_MC_Model, self).__init__()
        self.massconRelation = massconRelation

    def guess(self):

        guess = [10**(np.random.uniform(13, 16))]

        guess[0] = guess[0] / self.massScale

        return guess


    def paramLimits(self):

        return {'m200' : (1e13/self.massScale,1e16/self.massScale)}


    def __call__(self, x, m200):

        c200 = self.massconRelation(m200*self.massScale, self.zcluster, self.overdensity)        

        return super(NFW_MC_Model, self).__call__(x, m200, c200)

###############################

class NFWFitter(object):

    def __init__(self, profileBuilder, massconRelation):

        self.profileBuilder = profileBuilder

        if massconRelation is None:
            self.model = NFW_Model()
        else:
            self.model = NFW_MC_Model(massconRelation)


    ######

    def betaMethod(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster, guess = []):


        if guess == []:
            guess = self.model.guess()

        self.model.setData(beta_s, beta_s2, zcluster)

        fitter = fitmodel.FitModel(r_mpc, ghat, sigma_ghat, self.model,
                                   guess = guess)
        fitter.m.limits = self.model.paramLimits()
        fitter.fit()
        if fitter.have_fit:
            return fitter.par_vals
        return None

    #######

    def bootstrapFit(self, catalog, config):

        fitresults = []
        nfail = 0


        for i in range(config.nbootstraps):

            if i == 0:
                curCatalog = catalog
            else:
                curBootstrap = np.random.randint(0, len(catalog), size=len(catalog))
                curCatalog = catalog.filter(curBootstrap)

            beta_s = np.mean(curCatalog['beta_s'])
            beta_s2 = np.mean(curCatalog['beta_s']**2)

            r_mpc, ghat, sigma_ghat = self.profileBuilder(curCatalog, config)

            clean = sigma_ghat > 0

            fitresult = self.betaMethod(r_mpc[clean], ghat[clean], sigma_ghat[clean],  beta_s, beta_s2, 
                                   catalog.hdu.header['ZLENS'])

            if fitresult is None or not np.isfinite(m200):
                nfail += 1
            else:
                fitresults.append(fitresult)
                



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

    bootstrap_vals = fitter.bootstrapFit(catalog, config)

    savefit(bootstrap_vals, outputname)

############################

class FailedCreationException(Exception): pass


if __name__ == '__main__':


    catname = sys.argv[1]
    configname = sys.argv[2]
    outname = sys.argv[3]
    

    runNFWFit(catname, configname, outname)

    if not os.path.exists(outname):
        raise FailedCreationException


