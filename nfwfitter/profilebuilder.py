''' Set of utilities that read in simulation files and output shear profiles'''
##########################

import numpy as np

import nfwutils

#########################

def verifyShear(cat, col = 'gamma1_inf', thresh = 1.):

    assert(np.max(np.abs(cat[col])) < thresh)

#####


class ProfileBuilder(object):

    def configure(self, config):

        self.rescalecluster = config['rescalecluster']
        self.galaxypicker = config['galaxypicker']
        self.betacalcer = config['betacalcer']
        self.shearnoiser = config['shearnoiser']
        self.centergenerator = config['centergenerator']
        self.binner = config['binner']
        self.binnoiser = config['binnoiser']


    def __call__(self, sim):

        verifyShear(sim, 'gamma1_inf')
        verifyShear(sim, 'gamma2_inf')

        rescaledsim = self.rescalecluster(sim)

        verifyShear(rescaledsim, 'gamma1_inf')
        verifyShear(rescaledsim, 'gamma2_inf')


        galaxies = self.galaxypicker(rescaledsim)

        verifyShear(galaxies, 'gamma1_inf')
        verifyShear(galaxies, 'gamma2_inf')


        galaxies3d = self.betacalcer(galaxies)

        verifyShear(galaxies3d, 'gamma1_inf')
        verifyShear(galaxies3d, 'gamma2_inf')

        
        betas = galaxies3d.beta_s
        kappa = betas*galaxies3d.kappa_inf
        galaxies3d.g1 = betas*galaxies3d.gamma1_inf/(1 - kappa)
        galaxies3d.g2 = betas*galaxies3d.gamma2_inf/(1 - kappa)

        verifyShear(galaxies3d, 'g1')
        verifyShear(galaxies3d, 'g2')

        
        noisygalaxies = self.shearnoiser(galaxies3d)

        verifyShear(noisygalaxies, 'g1')
        verifyShear(noisygalaxies, 'g2')


        centeroffsetx, centeroffsety = self.centergenerator(noisygalaxies)
        print 'Center Offset:', centeroffsetx, centeroffsety

        delta_x = noisygalaxies.x_arcmin - centeroffsetx
        delta_y = noisygalaxies.y_arcmin - centeroffsety
        r_arcmin = np.sqrt(delta_x**2 + delta_y**2)


        dL = nfwutils.global_cosmology.angulardist(noisygalaxies.zlens)    

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

        
        noisygalaxies.r_arcmin = r_arcmin
        noisygalaxies.r_mpc = r_mpc
        noisygalaxies.ghat = E
        noisygalaxies.gcross = B

        verifyShear(noisygalaxies, 'ghat')
        verifyShear(noisygalaxies, 'gcross')

        
        
        profile = self.binner(noisygalaxies)

        verifyShear(noisygalaxies, 'ghat')


        clean = profile.sigma_ghat > 0
        cleanprofile = profile.filter(clean)
        cleanprofile.zcluster = noisygalaxies.zcluster
        cleanprofile.zlens = noisygalaxies.zlens

        noisyprofile = self.binnoiser(cleanprofile)

        verifyShear(noisyprofile, 'ghat')

        return noisyprofile
        

    


