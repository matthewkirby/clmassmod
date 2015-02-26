############
# Provides basic binning of galaxies into radial bins
###########

import numpy as np

#################################

class Binner(object):

    def __call__(self, catalog, config):

        maskedCat = catalog.filter(catalog['mask'])
        radii, shear, shearerr, avebeta, avebeta2, number = self._makeProfile(maskedCat, config)
        numbermask = number != -1

        if 'shearprofileerr' in config and config.shearprofileerr == 'gaussianapprox':
            allradii, allshear, allshearerr, allavebeta, allavebeta2, allnumber = self._makeProfile(catalog, config)
            assert(len(radii) == len(allradii))
            scalederr = allshearerr[numbermask]*np.sqrt(allnumber[numbermask]/number[numbermask])
            return radii[numbermask], shear[numbermask], scalederr

        return radii, shear, shearerr, avebeta, avebeta2

##################################

class dumbequalbins(Binner):

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
        

    def _makeProfile(self, catalog, config):

        sorted_cat = catalog.filter(np.argsort(catalog[self.profileCol]))
        sorted_cat = sorted_cat.filter(np.logical_and(sorted_cat[self.profileCol] > self.minradii, 
                                                      sorted_cat[self.profileCol] < self.maxradii))
        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []
        for i in range(0, len(sorted_cat), self.ngals):
            maxtake = min(i+self.ngals, len(catalog))
            radii.append(np.mean(sorted_cat['r_mpc'][i:maxtake]))
            shear.append(np.mean(sorted_cat['ghat'][i:maxtake]))
            shearerr.append(np.std(sorted_cat['ghat'][i:maxtake])/np.sqrt(maxtake-i))
            avebeta.append(np.mean(sorted_cat['beta_s'][i:maxtake]))
            avebeta2.append(np.mean(sorted_cat['beta_s'][i:maxtake]**2))
            ngals.append(len(sorted_cat['r_mpc'][i:maxtake]))

        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2), np.array(ngals)

############################

def bootstrapmean(distro, nboot=1000):

    bootedmeans = np.zeros(nboot)
    for i in range(nboot):
        curboot = np.random.randint(0, len(distro), len(distro))
        bootedmeans[i] = np.mean(distro[curboot])

    return np.mean(bootedmeans), np.std(bootedmeans)

############################

class bootstrapequalbins(Binner):

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
        

    def _makeProfile(self, catalog, config):

        sorted_cat = catalog.filter(np.argsort(catalog[self.profileCol]))
        sorted_cat = sorted_cat.filter(np.logical_and(sorted_cat[self.profileCol] > self.minradii, 
                                                      sorted_cat[self.profileCol] < self.maxradii))
        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []
        for i in range(0, len(sorted_cat), self.ngals):
            maxtake = min(i+self.ngals, len(catalog))
            radii.append(np.mean(sorted_cat['r_mpc'][i:maxtake]))

            curmean, curerr = bootstrapmean(sorted_cat['ghat'][i:maxtake])
            shear.append(curmean)
            shearerr.append(curerr)
            avebeta.append(np.mean(sorted_cat['beta_s'][i:maxtake]))
            avebeta2.append(np.mean(sorted_cat['beta_s'][i:maxtake]**2))
            ngals.append(len(sorted_cat['ghat'][i:maxtake]))

        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2), np.array(ngals)
            

##############################

class bootstrapfixedbins(Binner):

    def __init__(self, config = None):
        self.ngals = 200
        self.maxradii = 3.
        self.minradii = 0.
        self.binspacing = 'linear'
        self.nbins = 12.
        self.profileCol = 'r_mpc'

        if config is not None:
            self.maxradii = config.profilemax
            self.minradii = config.profilemin
            self.binspacing = config.binspacing
            self.nbins = config.nbins
            self.profileCol = config.profilecol
        

    def _makeProfile(self, catalog, config):

        if self.binspacing == 'linear':
            binedges = np.linspace(self.minradii, self.maxradii, self.nbins+1)
        else:
            binedges = np.logspace(np.log10(self.minradii), np.log10(self.maxradii), self.nbins+1)

        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []
        for i in range(self.nbins):
            mintake = binedges[i]
            maxtake = binedges[i+1]
            selected = catalog.filter(np.logical_and(catalog[self.profileCol] >= mintake,
                                                     catalog[self.profileCol] < maxtake))

        
        

            if len(selected) < 2:
                radii.append(-1)
                shear.append(-1)
                shearerr.append(-1)
                avebeta.append(-1)
                avebeta2.append(-1)
                ngals.append(-1)
                continue

            

            radii.append(np.mean(selected['r_mpc']))

            curmean, curerr = bootstrapmean(selected['ghat'])
            shear.append(curmean)
            shearerr.append(curerr)
            avebeta.append(np.mean(selected['beta_s']))
            avebeta2.append(np.mean(selected['beta_s']**2))
            ngals.append(len(selected))

        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2), np.array(ngals)
      


####################

class gaussianfixedbins(Binner):

    def __init__(self, config = None):
        self.ngals = 200
        self.maxradii = 3.
        self.minradii = 0.
        self.binspacing = 'linear'
        self.nbins = 12.
        self.profileCol = 'r_mpc'

        if config is not None:
            self.maxradii = config.profilemax
            self.minradii = config.profilemin
            self.binspacing = config.binspacing
            self.nbins = config.nbins
            self.profileCol = config.profilecol
            self.shapenoise = config.shapenoise
        

    def _makeProfile(self, catalog, config):

        if self.binspacing == 'linear':
            binedges = np.linspace(self.minradii, self.maxradii, self.nbins+1)
        else:
            binedges = np.logspace(np.log10(self.minradii), np.log10(self.maxradii), self.nbins+1)

        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []
        for i in range(self.nbins):
            mintake = binedges[i]
            maxtake = binedges[i+1]
            selected = catalog.filter(np.logical_and(catalog[self.profileCol] >= mintake,
                                                     catalog[self.profileCol] < maxtake))

            ngal = len(selected)            

            if ngal == 0:
                continue

            radii.append(np.mean(selected['r_mpc']))
            shear.append(np.mean(selected['ghat']))
            shearerr.append(self.shapenoise / np.sqrt(ngal))
            avebeta.append(np.mean(selected['beta_s']))
            avebeta2.append(np.mean(selected['beta_s']**2))
            ngals.append(ngal)


        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2), np.array(ngals)
      
