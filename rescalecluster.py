'''
Take a simulation, and pretend that it is at a different redshift

There's a few keys here:
1.) Substitute redshift
2.) rescale angle (keep physical scales fixed)
3.) rescale shear signal
'''

import numpy as np

import nfwutils

###############

class NoRedshiftRescaling(object):

    def __call__(self, sim):
        sim.zlens = sim.zcluster
        return sim


########

class RedshiftRescaler(object):

    def configure(self, config):
        self.targetz = config['targetz']

    def __call__(self, sim):

        newcat = sim.copy()
        newcat.zcluster = sim.zcluster
        newcat.zlens = self.targetz


        refDl = nfwutils.global_cosmology.angulardist(sim.zcluster)
        newDl = nfwutils.global_cosmology.angulardist(newcat.zlens)
        
        #rescale angle
        newcat.x_arcmin = (newcat.x_mpc/newDl)*(180.*60/np.pi)
        newcat.y_arcmin = (newcat.y_mpc/newDl)*(180.*60/np.pi)

        #rescale shear, kappa. Both are proportional to Dl. Beta_inf also changes.
        old_beta_inf = nfwutils.global_cosmology.beta([1e6], sim.zcluster)
        new_beta_inf = nfwutils.global_cosmology.beta([1e6], newcat.zlens)
        ratio = (newDl/refDl)*(new_beta_inf/old_beta_inf)  #this is in effect rescaling sigma_crit
        newcat.gamma1_inf = newcat.gamma1_inf*ratio
        newcat.gamma2_inf = newcat.gamma2_inf*ratio
        newcat.kappa_inf = newcat.kappa_inf*ratio

        return newcat
