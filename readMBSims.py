#######################
# Reads in MXXL files for processing
#######################

import numpy as np
import cPickle
import nfwutils
import nfwmodeltools as tools

#######################



#################################

class MBSimReader(object):

    def __init__(self, *args, **keywords):

        pass

    #########

    def getCosmology(self):

        return nfwutils.Cosmology(omega_m = 0.3, omega_l = 0.7, h=0.7)

    #########

    def load(self, filebase):

        return MBSim(filebase)

    
###############################


class MBSim(object):


    #########

    def __init__(self, filebase):

        # returns set of angle and mpc distances from the true center (flattened)
        # for each position, provide reduced shear g1 and g2, as well as z and beta for each source

        with open(filebase, 'rb') as input:
            clusterinfo = cPickle.load(input)

        m200 = clusterinfo['lm200']
        c200 = clusterinfo['lc200']

        nobjs = 10000

        delta_mpc = np.row_stack([np.random.uniform(-3, 3, size=nobjs),
                                  np.random.uniform(-3, 3, size=nobjs)])

        zcluster = 0.5
        Da = nfwutils.global_cosmology.angulardist(zcluster)

        delta_arcmin = (delta_mpc/Da)*(180./np.pi)*60.

        r_arcmin = np.sqrt(delta_arcmin[0]**2 + delta_arcmin[1]**2)
        r_mpc = np.sqrt(delta_mpc[0]**2 + delta_mpc[1]**2)

        cosphi = delta_mpc[0] / r_mpc
        sinphi = delta_mpc[1] / r_mpc
    
        sin2phi = 2.0*sinphi*cosphi
        cos2phi = 2.0*cosphi*cosphi-1.0

        self.zcluster = zcluster

        self.x_arcmin = delta_arcmin[0]
        self.y_arcmin = delta_arcmin[1]
        self.r_mpc = r_mpc
        self.r_arcmin = r_arcmin
        self.cos2phi = cos2phi
        self.sin2phi = sin2phi


        background_z = 2.
        self.redshifts = background_z*np.ones(nobjs)
        self.beta_s = nfwutils.global_cosmology.beta_s([background_z], self.zcluster)*np.ones_like(self.redshifts)

        
        r_scale = nfwutils.rscaleConstM(m200, c200, zcluster, 200)
        rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2
        
        gamma_t = tools.NFWShear(self.r_mpc, c200, r_scale, rho_c_over_sigma_c)
        kappa = tools.NFWShear(self.r_mpc, c200, r_scale, rho_c_over_sigma_c)

        g_t = self.beta_s*gamma_t/(1-self.beta_s*kappa)

        
        
        self.g1 = -g_t*self.cos2phi
        self.g2 = -g_t*self.sin2phi



        
        
    


    
    
    

    
    
    
