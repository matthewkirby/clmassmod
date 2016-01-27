''' Set of utilities that read in simulation files and output shear profiles'''
##########################

import numpy as np

import galaxypicker
import shearnoiser
import betacalcer
import centergenerator
import binnoiser

#########################

def BinnerFactory(config):

    binnermodule, binnerclass = config['binner'].split(':')
    binner = simutils.buildObject(binnermodule, binnerclass, config = config)

    return binnercalcer


########################

########################

def ProfileBuilderFactory(config):

    galaxypicker = galaxypicker.GalaxyPickerFactory(config) #could be multiple nested
    betacalcer = betacalcer.BetaCalcerFactory(config)
    shearnoiser = shearnoiser.ShearNoiserFactory(config) #could be multiple nested
    centergenerator = centergenerator.CenterGeneratorFactory(config)    
    binner = BinnerFactory(config)
    binnoiser = binnoiser.BinNoiserFactory(config)  #could be multiple nested

    

    return ProfileBuilder(galaxypicker, betacalcer, 
                          shearnoiser, centergenerator, 
                          binner, binnoiser)

###

class ProfileBuilder(object):

    def __init__(self, galaxypicker, betacalcer, 
                 shearnoiser, centergenerator, 
                 binner, binnoiser):

        self.galaxypicker = galaxypicker
        self.betacalcer = betacalcer
        self.shearnoiser = shearnoiser
        self.centergenerator = centergenerator
        self.binner = binner
        self.binnoiser = binnoiser


    def __call__(self, sim):

        galaxies = self.galaxypicker(sim)

        galaxies3d = self.betacalcer(galaxies)
        
        betas = galaxies3d.beta_s
        kappa = betas*galaxies3d.kappa_inf
        galaxies3d.g1 = betas*galaxies3d.gamma1_inf/(1 - kappa)
        galaxies3d.g2 = betas*galaxies3d.gamma2_inf/(1 - kappa)
        
        noisygalaxies = self.shearnoiser(galaxies3d)

        centeroffsetx, centeroffsety = self.centergenerator(sim)

        delta_x = noisygalaxies.x_arcmin - centeroffsetx
        delta_y = noisygalaxies.y_arcmin - centeroffsety
        r_arcmin = np.sqrt(delta_x**2 + delta_y**2)


        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    

        deltax_mpc = (delta_x * dL * np.pi)/(180.*60)
        deltay_mpc = (delta_y * dL * np.pi)/(180.*60)
        r_mpc = np.sqrt(deltax_mpc**2 + deltay_mpc**2)

        cosphi = deltax_mpc / r_mpc
        sinphi = deltay_mpc / r_mpc
        sin2phi = 2.0*sinphi*cosphi
        cos2phi = 2.0*cosphi*cosphi-1.0

        e1 = noisygalaxies.g1
        e2 = noisygalaxies.g2

        E = -(e1*cos2phi+e2*sin2phi)

        b1 =  e2
        b2 = -e1
        B = -(b1*cos2phi+b2*sin2phi)




        cols = [pyfits.Column(name = 'r_arcmin', format = 'E', array = r_arcmin),
                pyfits.Column(name = 'r_mpc', format='E', array = r_mpc),
                pyfits.Column(name = 'ghat', format='E', array = E),
                pyfits.Column(name = 'gcross', format='E', array = B),
                pyfits.Column(name = 'z', format='E', array = noisygalaxies.redshifts),
                pyfits.Column(name = 'beta_s', format = 'E', array = noisygalaxies.beta_s)]
        catalog = ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))
        catalog.hdu.header['ZLENS'] = sim.zcluster

        
        r_mpc, ghat, sigma_ghat, beta_s, beta_s2 = [x.astype(np.float64) for x in self.binner(catalog)]
        clean = sigma_ghat > 0


        profile = dict(r_mpc = r_mpc[clean], 
                       ghat = ghat[clean], 
                       sigma_ghat = sigma_ghat[clean], 
                       beta_s = beta_s[clean], 
                       beta_s2 = beta_s2[clean], 
                       zlens = sim.zcluster)

        noisyprofile = self.binnoiser(profile)

        return noisyprofile
        

    


