############
# Provides basic binning of galaxies into radial bins
###########

import numpy as np

class dumbequalbins(object):

    def __init__(self, config = None):
        self.ngals = 200
        self.maxradii = 2000
        self.minradii = 0
        self.profileCol = 'r_mpc'

        if config is not None:
            self.ngals = config.ngals
            self.maxradii = config.profilemax
            self.minradii = config.profilemin
            self.profileCol = config.profilecol
        

    def __call__(self, catalog, config):

        sorted_cat = catalog.filter(np.argsort(catalog[self.profileCol]))
        sorted_cat = sorted_cat.filter(np.logical_and(sorted_cat[self.profileCol] > self.minradii, 
                                                      sorted_cat[self.profileCol] < self.maxradii))
        radii = []
        shear = []
        shearerr = []
        for i in range(0, len(sorted_cat), self.ngals):
            maxtake = min(i+self.ngals, len(catalog))
            radii.append(np.mean(sorted_cat['r_mpc'][i:maxtake]))
            shear.append(np.mean(sorted_cat['ghat'][i:maxtake]))
            shearerr.append(np.std(sorted_cat['ghat'][i:maxtake])/np.sqrt(maxtake-i))

        return np.array(radii), np.array(shear), np.array(shearerr)
            

