'''
A facility to alter the binned shear profile after it has been created
'''

######

import numpy as np

import shearnoiser
import readtxtfile

######

class NoBinNoise(object):

    def __call__(self, profile):

        return profile
#####


        




