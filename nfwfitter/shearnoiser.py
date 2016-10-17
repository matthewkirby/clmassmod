''' Add noise & offsets to invidivual galaxy shear measurements'''
#############

import numpy as np

#############


class ShearNoiser(object):

    def __call__(self, sim):

        newg1, newg2 = self.addNoise(sim)

        newcat = sim.copy()

        newcat.g1 = newg1
        newcat.g2 = newg2

        return newcat

###############

class NoNoise(object):

    def __call__(self, sim):

        return sim

###############

class GaussianShapeNoise(ShearNoiser):

    def configure(self, config):

        self.shapenoise = config['shapenoise']

    def addNoise(self, sim):

        ngals = len(sim)
        newg1 = sim.g1 + self.shapenoise*np.random.standard_normal(ngals)
        newg2 = sim.g2 + self.shapenoise*np.random.standard_normal(ngals)

        return newg1, newg2
        

    
