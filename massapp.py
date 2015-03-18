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

    zcluster = catalog.hdu.header['ZLENS']
    dL = nfwutils.global_cosmology.angulardist(zcluster)

    r2 = (config.massappr2/3600.)*(np.pi/180.)*dL
    rmax = (config.massapprmax/3600.)*(np.pi/180.)*dL

    minradii = config.profilemin
    maxradii = config.profilemax
    nbins = config.nbins
    shapenoise = config.shapenoise

    # fit NFW profile

    nfwfitter = nfwfit.buildFitter(config)
    nfwm200, nfwm200err = nfwfitter.runUntilNotFail(catalog, config)
    nfwm200 = nfwm200['m200']
    c200 = nfwfitter.model.massconRelation(np.abs(nfwm200)*nfwfitter.model.massScale*nfwutils.global_cosmology.h, nfwfitter.model.zcluster, nfwfitter.model.overdensity)       


    rho_c = nfwutils.rho_crit(zcluster)
    rho_c_over_sigma_c = 1.5 * dL * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2

    nfwrscale = tools.rscaleConstM(nfwm200,
                                   c200,
                                   rho_c,
                                   200)
    

    # calculate gamma for catalog
    #use kappa from best fit nfw profile



    nfwkappa = tools.NFWKappa(catalog['r_mpc'], c200, nfwrscale, rho_c_over_sigma_c)
    gamma = catalog['ghat']*(1-catalog['beta']*nfwkappa)/catalog['beta']


    r1s = np.arange(minradii, maxradii, 0.1)
    kappa_proj = np.zeros_like(r1s)
    matching_m200s = np.zeros_like(r1s)
    mass_enclosed = np.zeros_like(r1s)
    density_enclosed = np.zeros_like(r1s)

    for cur_ap_index, r1 in enumerate(r1s):

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
        kappa_proj[cur_ap_index] = kappa_proj

        ##

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

        ##

        best_m200 = scipy.optimize.brentq(findNFW, 5e13, 1e16)
        matching_m200s[cur_ap_index] = best_m200
        best_c200 = nfwfitter.model.massconRelation(np.abs(best_m200)*nfwfitter.model.massScale*nfwutils.global_cosmology.h, nfwfitter.model.zcluster, nfwfitter.model.overdensity)       
        
        best_nfwrscale = tools.rscaleConstM(best_m200,
                                            best_c200,
                                            rho_c,
                                            200)

        mass_enclosed[cur_ap_index] = nfwutils.massInsideR(best_nfwrscale, best_c200,
                                                           zcluster, r1)
        vol = (4./3)*np.pi*r1**3
        density_enclosed[cur_ap_index] = mass_enclosed[cur_ap_index] / vol



    return r1s, kappa_proj, matching_m200s, mass_enclosed, density_enclosed


    
