import readtxtfile
import nfwutils
import numpy as np
import os.path

#####################

#__lss_dir__ = '/vol/euclid1/euclid1_raid1/schrabba/proj/spt/reduce201311/reduce_output_v3_ctim_bgq/massana/gb1_pofz_tim1_1p5m/{cluster}/massana_apera_c_x0/lssmocks'
#
#
#def loadLSSRealization(cluster, id = 'random'):
#        
#    clustermockdir = __lss_dir__.format(cluster = cluster)
#
#    if id == 'random':
#
#        lss_possible = readtxtfile.readtxtfile('{}/fit_lss.results'.format(clustermockdir))
#        ids_available = lss_possible[:,0][lss_possible[:,7] > 1000]
#
#        id = int(ids_available[np.random.randint(0, len(ids_available), 1)])
#
#    print 'Loading LSS %d' % id
#
#    lssfile = '{mockdir}/mock{id}/shear.profile.all.beta2.unchanged'.format(mockdir = clustermockdir,
#                                                                            id = id)
#
#    rawlssprofile = readtxtfile.readtxtfile(lssfile)
#
#    lssprofile = dict(r_mpc = rawlssprofile[:,0],
#                      gt = rawlssprofile[:,1],
#                      magbin = rawlssprofile[:,7])
#                      
#
#    return lssprofile
#
#
#
######################


class hstnoisebins(BinNoiser):
    '''This is a combination of a binner and a binnoiser. Do not combine with a shearnoiser.'''


    def __init__(self, config):

        super(hstnoisebins, self).__init__(config, *args, **kwds)

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
        self.betas = [x[5] for x in profile]
        self.beta2s = [x[6] for x in profile]
        self.magbinids = [x[-1] for x in profile]

        self.nbins = len(self.bincenters)

        self.useAveForCenter = False
        if 'centerforbin' in config and config['centerforbin'] == 'ave':
            self.useAveForCenter = True

        self.lssnoise = None
        if 'lssnoise' in config and config.lssnoise != 'False':
            self.lssnoise = config.lssnoise
        
            

    def __call__(self, catalog, config):

        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []

        lssrealization = None
        if self.lssnoise is not None:
            lssrealization = loadLSSRealization(self.clustername, self.lssnoise)


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
                lss_gt = lssrealization['gt'][selectbin]
                assert(lss_gt.shape == (1,))
                ghat += float(lss_gt)

            shear.append(ghat)  
        
            shearerr.append(self.deltag[i])
            avebeta.append(np.mean(selected['beta_s']))
            avebeta2.append(np.mean(selected['beta_s']**2))


        return np.array(radii), np.array(shear), np.array(shearerr), np.array(avebeta), np.array(avebeta2)
      
