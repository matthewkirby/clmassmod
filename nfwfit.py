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

def applyMask(x_arcmin, y_arcmin, config):

    def squaremask(x=0,y=0,theta=0, sidelength = 1.):
        #x,y in arcminutes

        theta_rad = np.pi*theta/180.
        
        dX = x_arcmin - x   
        dY = y_arcmin - y
        
        rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
        rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

        return (np.abs(rdX) + np.abs(rdY)) <= (np.sqrt(2)*sidelength/2.)

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
    cos2phi = sim.cos2phi
    sin2phi = sim.sin2phi

    E = -(e1*cos2phi+e2*sin2phi)

    b1 =  e2
    b2 = -e1
    B = -(b1*cos2phi+b2*sin2phi)

    if 'shapenoise' in config:
        E = E + config.shapenoise*np.random.standard_normal(len(E))
        B = B + config.shapenoise*np.random.standard_normal(len(B))


    visiblemask = np.ones(len(E)) == 1.
    if 'maskname' in config:
        visiblemask = applyMask(sim.x_arcmin, sim.y_arcmin, config)

    densitymask = np.ones(len(E)) == 1.
    if 'nperarcmin' in config:
        densitymask = applyDensityMask(sim.x_arcmin, sim.y_arcmin, sim.zcluster, config)

    mask = np.logical_and(visiblemask, densitymask)

    r_arcmin = sim.r_arcmin
    r_mpc = sim.r_mpc
    redshifts = sim.redshifts
    beta_s = sim.beta_s

    cols = [pyfits.Column(name = 'r_arcmin', format = 'E', array = r_arcmin[mask]),
            pyfits.Column(name = 'r_mpc', format='E', array = r_mpc[mask]),
            pyfits.Column(name = 'ghat', format='E', array = E[mask]),
            pyfits.Column(name = 'gcross', format='E', array = B[mask]),
            pyfits.Column(name = 'z', format='E', array = redshifts[mask]),
            pyfits.Column(name = 'beta_s', format = 'E', array = beta_s[mask])]
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

   return buildObject(config.readermodule, config.readerclass, config = config)

    
        
########################

class NFW_Model(object):

    def __init__(self):

        self.massScale = 1e14
        self.overdensity = 200

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

    def prepData(self, curCatalog, config):
        
        beta_s = np.mean(curCatalog['beta_s'])
        beta_s2 = np.mean(curCatalog['beta_s']**2)
        
        r_mpc, ghat, sigma_ghat = self.profileBuilder(curCatalog, config)

        clean = sigma_ghat > 0

        zlens = curCatalog.hdu.header['ZLENS']

        return r_mpc[clean], ghat[clean], sigma_ghat[clean], beta_s, beta_s2, zlens


    ######

    def betaMethod(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster, 
                   guess = [],
                   useSimplex=False):


        if guess == []:
            guess = self.model.guess()

        self.model.setData(beta_s, beta_s2, zcluster)

        fitter = fitmodel.FitModel(r_mpc, ghat, sigma_ghat, self.model,
                                   guess = guess)
        fitter.m.limits = self.model.paramLimits()
        fitter.fit(useSimplex = useSimplex)
        if fitter.have_fit:
            
            return fitter.par_vals
        return None


    #######

    def runUntilNotFail(self, catalog, config):

        r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = self.prepData(catalog, config)

        for i in range(config.nbootstraps):

            try:
                fitresult = self.betaMethod(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens)


                if fitresult is not None:
                    return fitresult

            except ValueError:
                pass

        #one last try with the SIMPLEX algorithm
        return self.betaMethod(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens, useSimplex=True)


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
                fitresult = self.betaMethod(r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens)

                

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

    fitvals = fitter.runUntilNotFail(catalog, config)
    
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


