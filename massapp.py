#!/usr/bin/env python
########
# My quick implementation of Henk's mass aperture scheme.
#########


import numpy as np
import nfwfit
import nfwutils
import nfwmodeltools as tools
import scipy.integrate
import scipy.optimize

##########



def avekappa(r1, r2, rscale, concentration, rho_c_over_sigma_c):

    def integrand(r):

        return r*tools.NFWKappa(r, concentration, rscale, rho_c_over_sigma_c)

    
    y,yerr = scipy.integrate.quad(integrand, r1, r2)
    avekappa = 2*y/(r2**2 - r1**2)

    return avekappa


#########

def logbinning(minradii, maxraii, nbins):

    binedges = np.logspace(np.log10(minradii), np.log10(maxradii), nbins+1)

    radii = []
    shear = []
    shearerr = []
    avebeta = []
    avebeta2 = []
    ngals = []
    for i in range(self.nbins):
        mintake = binedges[i]
        maxtake = binedges[i+1]
        selected = np.logical_and(catalog['r_mpc'] >= mintake,
                                  catalog['r_mpc'] < maxtake)

        ngal = len(gamma[selected])            

        if ngal == 0:
            continue

        radii.append(np.mean(catalog['r_mpc'][selected]))
        shear.append(np.mean(gamma[selected]))
        shearerr.append(self.shapenoise / np.sqrt(ngal))
        avebeta.append(np.mean(catalog['beta_s'][selected]))
        avebeta2.append(np.mean(catalog['beta_s'][selected]**2))
        ngals.append(ngal)
        
    radii = np.array(radii)
    shear = np.array(shear)
    shearerr = np.array(shearerr)
    avebeta = np.array(avebeta)
    avebeta2 = np.array(avebeta2)
    ngals  = np.array(ngals)

    return radii, shear, shearerr, avebeta, avebeta2, ngals




#########


def massapp(catalog, config):

    minradii = config.profilemin
    maxradii = config.profilemax
    nbins = config.nbins
    shapenoise = config.shapenoise

    # fit NFW profile

    nfwfitter = nfwfit.buildFitter(config)
    nfwm200, nfwm200err = nfwfitter.runUntilNotFail(catalog, config)
    c200 = nfwfitter.model.massconRelation(np.abs(m200)*nfwfitter.model.massScale*nfwutils.global_cosmology.h, nfwfitter.model.zcluster, nfwfitter.model.overdensity)       

    zcluster = catalog['ZLENS']

    rho_c = nfwutils.rho_crit(zcluster)
    rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2

    nfwrscale = tools.rscaleConstM(nfwm200,
                                   c200,
                                   rho_c,
                                   200)
    

    # calculate gamma for catalog
    #use kappa from best fit nfw profile



    nfwkappa = tools.NFWKappa(catalog['r_mpc'], c200, nfwrscale, rho_c_over_sigma_c)
    gamma = catalog['ghat']*(1-catalog['beta']*nfwkappa)/catalog['beta']


    #gamma integrals

    radii, shear, shearerr, avebeta, avebeta2, ngals = logbinning(r1, r2, int1bins)

    integrand1 = shear/radii
    res = scipy.integrate.simps(integrand1, radii)
    int1 = 2*res

    radii, shear, shearerr, avebeta, avebeta2, ngals = logbinning(r2, rmax, int2bins)
    integrand2 = shear/radii
    int2 = 2*rmax**2*scipy.integrate.simps(integrand2, radii)/(rmax**2 - r2**2)

    zeta_c = int1 + int2
    

    #kappa aperture

    kappa_ap = avekappa(r2, rmax, nfwrscale, c200, rho_c_over_sigma_c)


    #find best matched nfw that reproduces kappa core

    kappa_r1 = zeta_c + kappa_ap


    def findNFW(m200):

        c200 = nfwfitter.model.massconRelation(np.abs(m200)*nfwfitter.model.massScale*nfwutils.global_cosmology.h, nfwfitter.model.zcluster, nfwfitter.model.overdensity)       
        
        nfwrscale = tools.rscaleConstM(m200,
                                       c200,
                                       rho_c,
                                       200)

        avekappa = tools.aveEnclosedKappa(np.array([r1]),
                                          c200,
                                          nfwrscale,
                                          rho_c_over_sigma_c)
        return avekappa - kappa_r1

    m200_best = scipy.optimize.brentq(findNFW, 5e13, 1e16)


    return m200_best


    
