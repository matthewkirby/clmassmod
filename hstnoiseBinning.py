import readtxtfile
import nfwutils
import numpy as np

#####################

__lss_dir__ = '/vol/euclid1/euclid1_raid1/schrabba/proj/spt/reduce201311/reduce_output_v3_ctim_bgq/massana/gb1_pofz_tim1_1p5m/SPT-CLJ2359-5009/massana_apera_c_x0/lssmocks'
__lss_valid_ids__ = readtxtfile.readtxtfile('/vol/euclid1/euclid1_raid1/dapple/mxxlsims/shearprofiles/lss.valid_ids')[:,0]

def loadLSSRealization(id = None):

    if id is None or id == 'random':
        id = __lss_valid_ids__[np.random.randint(0, len(__lss_valid_ids__), 1)]

    lssfile = '%s/mock%d/shear.profile.all.beta2.unchanged' % (__lss_dir__,
                                                               id)

    rawlssprofile = readtxtfile.readtxtfile(lssfile)

    lssprofile = dict(r_mpc = rawlssprofile[:,0],
                      gt = rawlssprofile[:,1],
                      magbin = rawlssprofile[:,7])
                      

    return lssprofile



#####################


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
        self.magbinids = [x[-1] for x in profile]

        self.nbins = len(self.bincenters)

        self.useAveForCenter = False
        if 'centerforbin' in config and config['centerforbin'] == 'ave':
            self.useAveForCenter = True

        self.lssnoise = None
        if 'lssnoise' in config:
            self.lssnoise = config.lssnoise
        
            

    def __call__(self, catalog, config):

        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []

        lssrealization = None
        if self.lssnoise:
            lssrealization = loadLSSRealization(self.lssnoise)


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
                
            #Take the mean shear and add noise
            ghat = np.mean(selected['ghat']) + self.deltag[i]*np.random.standard_normal()

            #if applicable, add LSS noise
            if lssrealization is not None:
                selectbin = np.logical_and(lssrealization['r_mpc'] == self.bincenters[i],
                                           lssrealization['magbin'] == self.magbinids[i])
                ghat += lssrealization['gt'][selectbin]

            shear.append(ghat)  
        
            shearerr.append(self.deltag[i])
            avebeta.append(np.mean(selected['beta_s']))
            avebeta2.append(np.mean(selected['beta_s']**2))


        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2)
      
