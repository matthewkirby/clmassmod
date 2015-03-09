import readtxtfile
import nfwutils
import numpy as np

class hstnoisebins(object):
    ''' Note: This binning class adds noise, unlike the other binning classes. Should not be run with shape noise added at the catalog reader stage.'''

    def __init__(self, config):

        if 'shapenoise' in config and config['shapenoise'] > 0.:
            raise ValueError
        
        self.maxradii = config.profilemax
        self.minradii = config.profilemin

        self.profileCol = config.profilecol
        self.binwidth = config.binwidth
        self.profilefile = config.profilefile

        profile = readtxtfile.readtxtfile(self.profilefile)
        
        self.bincenters = [x[0] for x in profile]
        self.deltag = [x[2] for x in profile]

        self.nbins = len(self.bincenters)

        self.useAveForCenter = False
        if 'centerforbin' in config and config['centerforbin'] == 'ave':
            self.useAveForCenter = True
        
            

    def __call__(self, catalog, config):

        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []


        for i in range(self.nbins):

            if self.bincenters[i] < self.minradii or self.bincenters[i] > self.maxradii:
                continue

            selected = catalog.filter(np.logical_and(catalog[self.profileCol] >= (self.bincenters[i] - self.binwidth/2.),
                                                     catalog[self.profileCol] < (self.bincenters[i] + self.binwidth/2.)))

            ngal = len(selected)            

            if ngal == 0:
                continue

            

            if self.useAveForCenter:
                radii.append(np.mean(selected[self.profileCol]))
            else:
                radii.append(self.bincenters[i])


            shear.append(np.mean(selected['ghat']) + self.deltag[i]*np.random.standard_normal())  #Take the mean shear and add noise
        
            shearerr.append(self.deltag[i])
            avebeta.append(np.mean(selected['beta_s']))
            avebeta2.append(np.mean(selected['beta_s']**2))


        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2)
      
