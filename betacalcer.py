'''
Assign betas and redshifts to galaxies, and assign appropriately.
'''

##########

import simutils
import numpy as np
import nfwutils

##########


class BetaCalcer(object):
    '''ABSTRACT'''

    def configure(self, config):

        self.targetz = None
        if 'targetz' in config:
            self.targetz = config['targetz']
        

    def __call__(self, galaxies):


        zlens = galaxies.zcluster
        if self.targetz is not None:
            zlens = self.targetz

        beta_s = self.calcBetas(zlens, galaxies)

        if self.targetz is not None:
            Dl_ref = nfwutils.global_cosmology.angulardist(galaxies.zcluster)
            Dl_target = nfwutils.global_cosmology.angulardist(self.targetz)
            beta_s = beta_s*(Dl_target/Dl_ref)


        galaxies.beta_s = beta_s

        return galaxies

#########

class OnlyRescaleLens(BetaCalcer):

    def calcBetas(zlens, galaxies):
        return np.ones(len(galaxies))

#########


class FixedRedshift(BetaCalcer):

    def configure(self, config):

        super(FixedRedshift, self).configure(config)
        self.z_source = config['zsource']


    def calcBetas(self, zlens, galaxies):

        beta_s = nfwutils.global_cosmology.beta_s(np.ones(len(galaxies))*self.z_source, zlens)

        return beta_s


##########

class FixedBeta(BetaCalcer):

    def configure(self, config):

        super(FixedBeta, self).configure(config)
        self.beta_s = config['beta_s']

    def calcBetas(self, zlens, galaxies):


        return self.beta_s*np.ones(len(galaxies))

        
        
