''' Add noise & offsets to invidivual galaxy shear measurements'''
#############

import simutils


#############


class ShearNoiser(object):

    def __call__(self, sim):

        noisygals = self.prevnoiser(sim)

        newg1, newg2 = self.addNoise(noisygals)

        noisygals.g1 = newg1
        noisygals.g2 = newg2

        return noisygals

###############

class NoNoise(object):

    def __call__(self, sim):

        return sim

###############

class GaussianShapeNoise(ShearNoiser):

    def configure(self, config):

        self.shapenoise = config['shapenoise']

    def addNoise(sim):

        ngals = len(sim)
        newg1 = sim.g1 + self.shapenoise*np.random.standard_normal(ngals)
        newg2 = sim.g2 + self.shapenoise*np.random.standard_normal(ngals)

        return newg1, newg2
        

    
