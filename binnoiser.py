'''
A facility to alter the binned shear profile after it has been created
'''

######

import simutils

######


def BinNoiserFactory(config):
    '''Build a BinNoiser object from a config file. Could be multiple, nested.'''

    binnoiser = NoBinNoise()    
    if 'binnoisers' in config:
        for nextnoiser in config['binnoisers'].split(','):
            noisermodule, noiserclass = nextnoiser.split(':')
            binnoiser = simutils.buildObject(noisermodule, 
                                             noiserclass, 
                                             prevnoiser = binnoiser, 
                                             config = config)

    return picker

#####

class NoBinNoise(object):

    def __call__(self, profile):

        return profile
#####

class BinNoiser(object):

    def __init__(self, prevnoiser, config):

        self.prevnoiser = prevnoiser
        self.config = config

    def __call__(self, profile):

        prevnoiseprofile = self.prevnoiser(profile)
        
        return self.addNoise(prevnoiseprofile)




