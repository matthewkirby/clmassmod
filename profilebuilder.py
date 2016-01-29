''' Set of utilities that read in simulation files and output shear profiles'''
##########################

import numpy as np


#########################


class ProfileBuilder(object):

    def configure(self, config):

        self.galaxypicker = config['galaxypicker']
        self.betacalcer = config['betacalcer']
        self.shearnoiser = config['shearnoiser']
        self.centergenerator = config['centergenerator']
        self.binner = config['binner']
        self.binnoiser = config['binnoiser']


    def __call__(self, sim):

        galaxies = self.galaxypicker(sim)

        galaxies3d = self.betacalcer(galaxies)
        
        betas = galaxies3d.beta_s
        kappa = betas*galaxies3d.kappa_inf
        galaxies3d.g1 = betas*galaxies3d.gamma1_inf/(1 - kappa)
        galaxies3d.g2 = betas*galaxies3d.gamma2_inf/(1 - kappa)
        
        noisygalaxies = self.shearnoiser(galaxies3d)

        centeroffsetx, centeroffsety = self.centergenerator(noisygalaxies)

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

        
        noisygalaxies.r_arcmin = r_arcmin
        noisygalaxies.r_mpc = r_mpc
        noisygalaxies.ghat = E
        noisygalaxies.gcross = B
        
        
        profile = self.binner(noisygalaxies)
        clean = sigma_ghat > 0

        cleanprofile = profile.filter(clean)
        cleanprofile.zcluster = sim.zcluster
        noisyprofile = self.binnoiser(cleanprofile)

        return noisyprofile
        

    


