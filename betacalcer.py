'''
Assign betas and redshifts to galaxies, and assign appropriately.
'''

##########

import simutils
import numpy as np
import nfwutils

##########

def BetaCalcerFactory(config):



    betamodule, betaclass = config['betacalcer'].split(':')
    betacalcer = simutils.buildObject(betamodule, betaclass, config = config)

    return betacalcer

###########

class BetaCalcer(object):

    def __init__(self, config):
        self.config = config

    def __call__(self, galaxies):


        zlens = galaxies.zcluster
        if 'targetz' in self.config:
            zlens = self.config['targetz']

        beta_s = self.calcBetas(zlens, galaxies)

        if 'targetz' in self.config:
            Dl_ref = nfwutils.global_cosmology.angulardist(galaxies.zcluster)
            Dl_target = nfwutils.global_cosmology.angulardist(self.config['targetz'])
            beta_s = beta_s*(Dl_target/Dl_ref)


        galaxies.beta_s = beta_s

        return galaxies

#########

class OnlyRescaleLens(BetaCalcer):

    def calcBetas(zlens, galaxies):
        return np.ones(len(galaxies))

#########


class FixedRedshift(BetaCalcer):


    def calcBetas(self, zlens, galaxies):

        z_source = self.config.zsource

        beta_s = nfwutils.global_cosmology.beta_s(np.ones(len(galaxies))*z_source, zlens)

        return beta_s


##########

class FixedBeta(BetaCalcer):

    def calcBetas(self, zlens, galaxies):

        beta_s = self.config.beta_s
        return beta_s*np.ones(len(galaxies))

        
        
