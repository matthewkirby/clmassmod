''' Add noise & offsets to invidivual galaxy shear measurements'''
#############

import simutils


#############

def ShearNoiserFactory(config):
    '''Build a noiser object from a config file. Could be multiple, nested.'''

    noiser = NoNoise()    
    if 'shearnoisers' in config:
        for nextnoiser in config['shearnoisers'].split(','):
            noisermodule, noiserclass = nextpicker.split(':')
            noiser = simutils.buildObject(noisermodule, noiserclass, prevnoiser = noiser, config = config)

    return noiser

###############

class ShearNoiser(object):

    def __init__(self, prevnoiser, config):
        
        self.prevnoiser = prevnoiser
        self.config = config

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

    def addNoise(sim):

        shapenoise = self.config.shapenoise
        ngals = len(sim)
        newg1 = sim.g1 + shapenoise*np.random.standard_normal(ngals)
        newg2 = sim.g2 + shapenoise*np.random.standard_normal(ngals)

        return newg1, newg2
        

    
