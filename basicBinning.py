############
# Provides basic binning of galaxies into radial bins
###########

import numpy as np
import catalog

#################################

# Contents of a binnoiser

#        if 'shearprofileerr' in config and config.shearprofileerr == 'gaussianapprox':
#            allradii, allshear, allshearerr, allavebeta, allavebeta2, allnumber = self._makeProfile(catalog, config)
#            assert(len(radii) == len(allradii))
#            scalederr = allshearerr[numbermask]*np.sqrt(allnumber[numbermask]/number[numbermask])
#            return radii[numbermask], shear[numbermask], scalederr
#

#################################


class DumbEqualBins(object):

    def configure(self, config):

        self.ngals = 200
        self.maxradii = 2000
        self.minradii = 0
        self.profileCol = 'r_mpc'

        if 'ngals' in config:
            self.ngals = config['ngals']
            self.maxradii = config['profileMax']
            self.minradii = config['profileMin']
            self.profileCol = config['profilecol']
        

    def __call__(self, cat):

        profileCol = getattr(cat, self.profileCol)

        sorted_cat = cat.filter(np.argsort(profileCol))
        sorted_cat = sorted_cat.filter(np.logical_and(sorted_cat[self.profileCol] > self.minradii, 
                                                      sorted_cat[self.profileCol] < self.maxradii))
        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []
        for i in range(0, len(sorted_cat), self.ngals):
            maxtake = min(i+self.ngals, len(cat))
            radii.append(np.mean(sorted_cat[self.profileCol][i:maxtake]))
            shear.append(np.mean(sorted_cat['ghat'][i:maxtake]))
            shearerr.append(np.std(sorted_cat['ghat'][i:maxtake])/np.sqrt(maxtake-i))
            avebeta.append(np.mean(sorted_cat['beta_s'][i:maxtake]))
            avebeta2.append(np.mean(sorted_cat['beta_s'][i:maxtake]**2))
            ngals.append(len(sorted_cat[self.profileCol][i:maxtake]))

        profile = catalog.Catalog()
        setattr(profile, self.profileCol, np.array(radii))
        profile.ghat = np.array(shear)
        profile.sigma_ghat = np.array(shearerr)
        profile.beta_s = np.array(avebeta)
        profile.beta_s2 = np.array(avebeta2)
        profile.ngals = np.array(ngals)
        

        return profile

############################

def bootstrapmean(distro, nboot=1000):

    bootedmeans = np.zeros(nboot)
    for i in range(nboot):
        curboot = np.random.randint(0, len(distro), len(distro))
        bootedmeans[i] = np.mean(distro[curboot])

    return np.mean(bootedmeans), np.std(bootedmeans)

############################

class BootstrapEqualBins(object):

    def configure(self, config):
        self.ngals = 200
        self.maxradii = 2000
        self.minradii = 0
        self.profileCol = 'r_mpc'

        if 'ngals' in config:
            self.ngals = config['ngals']
            self.maxradii = config['profileMax']
            self.minradii = config['profileMin']
            self.profileCol = config['profilecol']
        

    def __call__(self, cat):

        profileCol = getattr(cat, self.profileCol)

        sorted_cat = cat.filter(np.argsort(profileCol))
        sorted_cat = sorted_cat.filter(np.logical_and(sorted_cat[self.profileCol] > self.minradii, 
                                                      sorted_cat[self.profileCol] < self.maxradii))
        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []
        for i in range(0, len(sorted_cat), self.ngals):
            maxtake = min(i+self.ngals, len(cat))
            radii.append(np.mean(sorted_cat[self.profileCol][i:maxtake]))

            curmean, curerr = bootstrapmean(sorted_cat['ghat'][i:maxtake])
            shear.append(curmean)
            shearerr.append(curerr)
            avebeta.append(np.mean(sorted_cat['beta_s'][i:maxtake]))
            avebeta2.append(np.mean(sorted_cat['beta_s'][i:maxtake]**2))
            ngals.append(len(sorted_cat['ghat'][i:maxtake]))

        profile = catalog.Catalog()
        setattr(profile, self.profileCol, np.array(radii))
        profile.ghat = np.array(shear)
        profile.sigma_ghat = np.array(shearerr)
        profile.beta_s = np.array(avebeta)
        profile.beta_s2 = np.array(avebeta2)
        profile.ngals = np.array(ngals)

        return profile

            

##############################

class BootstrapFixedBins(object):

    def configure(self, config):
        self.ngals = 200
        self.maxradii = 3.
        self.minradii = 0.
        self.binspacing = 'linear'
        self.nbins = 12
        self.profileCol = 'r_mpc'

        if 'nbins' in config:

            self.maxradii = config['profileMax']
            self.minradii = config['profileMin']
            self.binspacing = config['binspacing']
            self.nbins = config['nbins']
            self.profileCol = config['profilecol']
        

    def __call__(self, cat):

        profileCol = getattr(cat, self.profileCol)

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
            selected = cat.filter(np.logical_and(profileCol >= mintake,
                                                     profileCol < maxtake))

        
        

            if len(selected) < 2:
                radii.append(-1)
                shear.append(-1)
                shearerr.append(-1)
                avebeta.append(-1)
                avebeta2.append(-1)
                ngals.append(-1)
                continue

            

            radii.append(np.mean(getattr(selected,self.profileCol)))

            curmean, curerr = bootstrapmean(selected.ghat)
            shear.append(curmean)
            shearerr.append(curerr)
            avebeta.append(np.mean(selected.beta_s))
            avebeta2.append(np.mean(selected.beta_s**2))
            ngals.append(len(selected))

        profile = catalog.Catalog()
        setattr(profile, self.profileCol, np.array(radii))
        profile.ghat = np.array(shear)
        profile.sigma_ghat = np.array(shearerr)
        profile.beta_s = np.array(avebeta)
        profile.beta_s2 = np.array(avebeta2)
        profile.ngals = np.array(ngals)

        return profile

      


####################

class GaussianFixedBins(object):

    def configure(self, config):
        self.ngals = 200
        self.maxradii = 3.
        self.minradii = 0.
        self.binspacing = 'linear'
        self.nbins = 12.
        self.profileCol = 'r_mpc'

        if 'nbins' in config:
            self.maxradii = config['profileMax']
            self.minradii = config['profileMin']
            self.binspacing = config['binspacing']
            self.nbins = config['nbins']
            self.profileCol = config['profilecol']
            self.shapenoise = config['shapenoise']
        

    def __call__(self, cat):

        profileCol = getattr(cat, self.profileCol)

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
            selected = cat.filter(np.logical_and(profileCol >= mintake,
                                                     profileCol < maxtake))

            ngal = len(selected)            

            if ngal == 0:
                continue

            radii.append(np.mean(getattr(selected, self.profileCol)))
            shear.append(np.mean(selected.ghat))
            shearerr.append(self.shapenoise / np.sqrt(ngal))
            avebeta.append(np.mean(selected.beta_s))
            avebeta2.append(np.mean(selected.beta_s**2))
            ngals.append(ngal)



        profile = catalog.Catalog()
        setattr(profile, self.profileCol, np.array(radii))
        profile.ghat = np.array(shear)
        profile.sigma_ghat = np.array(shearerr)
        profile.beta_s = np.array(avebeta)
        profile.beta_s2 = np.array(avebeta2)
        profile.ngals = np.array(ngals)

        return profile






########################################


