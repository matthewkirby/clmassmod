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

    def __call__(self, galaxies):


        zlens = galaxies.zlens

        beta_s = self.calcBetas(zlens, galaxies)

        galaxies.beta_s = beta_s

        return galaxies

#########

class InfiniteRedshift(BetaCalcer):

    def calcBetas(self, zlens, galaxies):

        return np.ones(len(galaxies))

#########


class FixedRedshift(BetaCalcer):

    def configure(self, config):

        self.z_source = config['zsource']


    def calcBetas(self, zlens, galaxies):

        beta_s = nfwutils.global_cosmology.beta_s(np.ones(len(galaxies))*self.z_source, zlens)

        return beta_s


##########

class FixedBeta(BetaCalcer):

    def configure(self, config):

        self.beta = config['beta']

    def calcBetas(self, zlens, galaxies):


        return self.beta*np.ones(len(galaxies))/nfwutils.global_cosmology.beta([1e6], galaxies.zlens)

        
        
