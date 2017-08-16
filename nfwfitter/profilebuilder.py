''' Set of utilities that read in simulation files and output shear profiles'''
##########################

import numpy as np

import nfwutils

#########################

def verifyShear(cat, col = 'gamma1_inf', thresh = 5.):

    assert(np.max(np.abs(cat.table[col])) < thresh)

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

        rescaledsim = self.rescalecluster(sim)

        galaxies = self.galaxypicker(rescaledsim)


        galaxies3d = self.betacalcer(galaxies)



        
        betas = galaxies3d.beta_s
        kappa = betas*galaxies3d.kappa_inf
        galaxies3d.g1 = betas*galaxies3d.gamma1_inf/(1 - kappa)
        galaxies3d.g2 = betas*galaxies3d.gamma2_inf/(1 - kappa)
        g = np.sqrt(galaxies3d.g1**2 + galaxies3d.g2**2)
        no_arcs = np.abs(g) < 5
        galaxies_noarcs = galaxies3d.filter(no_arcs)

        verifyShear(galaxies_noarcs, 'g1')
        verifyShear(galaxies_noarcs, 'g2')

        
        noisygalaxies = self.shearnoiser(galaxies_noarcs)


        centeroffsetx, centeroffsety = self.centergenerator(noisygalaxies)
        print 'Center Offset:', centeroffsetx, centeroffsety

        delta_x = noisygalaxies.x_arcmin - centeroffsetx
        delta_y = noisygalaxies.y_arcmin - centeroffsety
        r_arcmin = np.sqrt(delta_x**2 + delta_y**2)


        dL = nfwutils.global_cosmology.angulardist(noisygalaxies.zlens)    

        deltax_mpc = (delta_x * dL * np.pi)/(180.*60)
        deltay_mpc = (delta_y * dL * np.pi)/(180.*60)
        r_mpc = np.sqrt(deltax_mpc**2 + deltay_mpc**2)

        # Convert gamma1 and gamma2 to tangential and cross shear terms
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




        clean = profile.sigma_ghat > 0
        cleanprofile = profile.filter(clean)
        cleanprofile.zcluster = noisygalaxies.zcluster
        cleanprofile.zlens = noisygalaxies.zlens

        noisyprofile = self.binnoiser(cleanprofile)



        return noisyprofile
        

    


