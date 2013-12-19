#!/usr/bin/env python
#######################
# Runs a bootstrapped NFW model fit to simulated shear data.
# Galaxies are binned into an average shear profile before fitting to NFW model.
# Options to explore radial fit range, mass-concentration relation, and binning scheme in fits.
########################

import importlib, cPickle, sys, os
import numpy as np
import pyfits
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

def readSimCatalog(catalogname, config=None):

    sim = ldac.LDACCat(pyfits.open(catalogname)[1])

    boxlength = sim['BOXWIDTHCOMOVINGHINVMPC'] / nfwutils.global_cosmology.h
    gridsize = np.sqrt(float(sim['A00'].shape[1]))

    clusterz = float(sim['ZLENS'])

    ####becker angles

    Ng = 512

    dc_readout = nfwutils.global_cosmology.comovingdist(clusterz)+sim['BOXLENGTHCOMOVINGHINVMPC'][0]/2.0
    da = np.arctan2(sim['BOXWIDTHCOMOVINGHINVMPC'][0],dc_readout)/np.pi*180.0*60.0/Ng

    xcen = Ng/2.0
    ycen = Ng/2.0

    radii_arcmin = np.zeros((Ng, Ng))
    cosphi = np.zeros((Ng,Ng))
    sinphi = np.zeros((Ng,Ng))

    for i in range(Ng):
        for j in range(Ng):
            xr = (i + 0.5 - xcen)*da
            yr = (j + 0.5 - ycen)*da
            r = np.sqrt(xr*xr + yr*yr)
            radii_arcmin[i,j] = r

            cosphi[i,j] = xr/r
            sinphi[i,j] = yr/r
    

    A00 = sim['A00'].reshape(gridsize, -1)
    A11 = sim['A11'].reshape(gridsize, -1)
    A10 = sim['A10'].reshape(gridsize, -1)
    A01 = sim['A01'].reshape(gridsize, -1)

    kappa = 0.5*(2- A00 - A11)
    gamma1 = 0.5*(A11 - A00)
    gamma2 = -0.5*(A01+A10)
    reducedShear1 = gamma1/(1-kappa)
    reducedShear2 = gamma2/(1-kappa)

    r_arcmin = radii_arcmin.flatten()

    r_mpc = (r_arcmin/60.)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(clusterz)

    cosphi = cosphi.flatten()
    sinphi = sinphi.flatten()
    sin2phi = 2.0*sinphi*cosphi
    cos2phi = 2.0*cosphi*cosphi-1.0
    

    e1 = reducedShear1.flatten()
    e2 = reducedShear2.flatten()

    E = -(e1*cos2phi+e2*sin2phi)

    b1 =  e2
    b2 = -e1
    B = -(b1*cos2phi+b2*sin2phi)

    redshifts = np.ones(len(r_mpc))*float(sim['ZSOURCE'])
    beta_s = nfwutils.global_cosmology.beta_s(redshifts, clusterz)

    cols = [pyfits.Column(name = 'r_arcmin', format = 'E', array = r_arcmin),
            pyfits.Column(name = 'r_mpc', format='E', array = r_mpc),
            pyfits.Column(name = 'ghat', format='E', array = E),
            pyfits.Column(name = 'gcross', format='E', array = B),
            pyfits.Column(name = 'z', format='E', array = redshifts),
            pyfits.Column(name = 'beta_s', format = 'E', array = beta_s)]
    catalog = ldac.LDACCat(pyfits.new_table(pyfits.ColDefs(cols)))
    catalog.hdu.header.update('ZLENS', clusterz)
    catalog.hdu.header.update('M200C', float(sim['M200C']))
    catalog.hdu.header.update('C200C', float(sim['C200C']))

            
    


    return catalog


########################

def readConfiguration(configname):

    config = bashreader.parseFile(configname)
    return config
    

########################

def buildObject(modulename, classname, *args, **kwds):

    print args,kwds

    aModule = importlib.import_module(modulename)
    aClass = getattr(aModule, classname)
    anObject = aClass(*args, **kwds)

    return anObject
    
#####

def buildFitter(config):

    profileBuilder = buildObject(config.profilemodule, config.profilebuilder, config)

    massconRelation = buildObject(config.massconmodule, config.massconrelation, config)

    fitter = NFWFitter(profileBuilder = profileBuilder, massconRelation = massconRelation)


    return fitter
        
########################


#####

class NFWFitter(object):

    def __init__(self, profileBuilder, massconRelation):

        self.profileBuilder = profileBuilder
        self.massconRelation = massconRelation


    ######

    def betaMethod(self, r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zcluster, guess = []):

        ####

        class NFWModel:
            def __init__(self, beta_s, beta_s2, zcluster, massconRelation, massScale):

                self.beta_s = beta_s
                self.beta_s2 = beta_s2
                self.zcluster = zcluster
                self.massconRelation = massconRelation
                self.massScale = massScale
                self.overdensity = 200


            def __call__(self, x, m200):
                
                concentration = self.massconRelation(m200*self.massScale, self.zcluster, self.overdensity)

                r_scale = nfwutils.rscaleConstM(m200*self.massScale, concentration, self.zcluster, self.overdensity)
                
                nfw_shear_inf = tools.NFWShear(x, concentration, r_scale, self.zcluster)
                nfw_kappa_inf = tools.NFWKappa(x, concentration, r_scale, self.zcluster)
                
                g = self.beta_s*nfw_shear_inf / (1 - ((self.beta_s2/self.beta_s)*nfw_kappa_inf) )
                
                return g

        ###

        massScale = 1e14

        if guess == []:
            guess = [10**(np.random.uniform(13, 16))]

        guess[0] = guess[0] / massScale

        model = NFWModel(beta_s, beta_s2, zcluster, self.massconRelation, massScale)

        fitter = fitmodel.FitModel(r_mpc, ghat, sigma_ghat, model,
                                   guess = guess)
        fitter.m.limits = {'m200' : (1e13/massScale,1e16/massScale)}
        fitter.fit()
        m200 = None
        if fitter.have_fit:
                       
            m200 = float(fitter.par_vals['m200'])*massScale

        return m200
        





        


    #######

    def bootstrapFit(self, catalog, config):

        mass = []
        nfail = 0
        for i in range(config.nbootstraps):


            curBootstrap = np.random.randint(0, len(catalog), size=len(catalog))

            curCatalog = catalog.filter(curBootstrap)

            beta_s = np.mean(curCatalog['beta_s'])
            beta_s2 = np.mean(curCatalog['beta_s']**2)

            r_mpc, ghat, sigma_ghat = self.profileBuilder(curCatalog, config)

            clean = sigma_ghat > 0

            m200 = self.betaMethod(r_mpc[clean], ghat[clean], sigma_ghat[clean],  beta_s, beta_s2, 
                                   catalog.hdu.header['ZLENS'])

            if m200 is None or not np.isfinite(m200):
                nfail += 1
            else:
                mass.append(m200)



        return np.array(mass), nfail




########################

def savefit(bootstrap_vals, outputname):

    with open(outputname, 'wb') as output:

        ### !!!!!!!!!!!!  Factor out h!
        cPickle.dump(np.array(bootstrap_vals)*nfwutils.global_cosmology.h, output, -1)


########################

def runNFWFit(catalogname, configname, outputname):


    config = readConfiguration(configname)

    catalog = readSimCatalog(catalogname, config)

    fitter = buildFitter(config)

    bootstrap_vals = fitter.bootstrapFit(catalog, config)

    savefit(bootstrap_vals, outputname)

############################

class FailedCreationException(Exception): pass


if __name__ == '__main__':

    cosmo = nfwutils.Cosmology(omega_m = 0.27, omega_l = 0.73, h = 0.7)
    nfwutils.global_cosmology.set_cosmology(cosmo)

    catname = sys.argv[1]
    configname = sys.argv[2]
    outname = sys.argv[3]
    

    runNFWFit(catname, configname, outname)

    if not os.path.exists(outname):
        raise FailedCreationException


