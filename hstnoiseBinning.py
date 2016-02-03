import numpy as np
import os.path


import readtxtfile
import nfwutils
import basicBinning
import betacalcer
import shearnoiser
import catalog

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

class HSTBinning(object):

    def configure(self, config):

        assert(isinstance(config['betacalcer'], betacalcer.InfiniteRedshift))
        assert(isinstance(config['shearnoiser'], shearnoiser.NoNoise))

        


        self.profilefile = config['profilefile']

        self.maxradii = config['profilemax']
        self.minradii = config['profilemin']
        self.binwidth = config['binwidth']

        self.profileCol = config['profilecol']



        profile = readtxtfile.readtxtfile(self.profilefile)
        
        self.bincenters = np.array([x[0] for x in profile])
        self.deltag = np.array([x[2] for x in profile])
        self.betas = np.array([x[5] for x in profile]) #not scaled by beta_inf

        self.magbinids = np.array([x[-1] for x in profile])

        self.nbins = len(self.bincenters)

        self.useAveForCenter = False
        if 'centerforbin' in config and config['centerforbin'] == 'ave':
            self.useAveForCenter = True


    ####

    def doBinning(self, galaxies):

        radii = []
        shear = []
        shearerr = []
        avebeta = []
        avebeta2 = []
        ngals = []

        profileCol = getattr(galaxies,self.profileCol)

        for i in range(self.nbins):

            mintake = self.bincenters[i] - self.binwidth/2.
            maxtake = self.bincenters[i] + self.binwidth/2.
            selected = galaxies.filter(np.logical_and(profileCol >= mintake,
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

            curmean, curerr = basicBinning.bootstrapmean(selected.ghat)
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


    #####


    def addNoise(self, profile):

        noisy = profile.copy()

        #rescale beta
        beta_s = self.betas/nfwutils.global_cosmology.beta([1e6], profile.zlens)
        newghat = profile.ghat*beta_s
        
        noisy.ghat = newghat + self.deltag*np.random.standard_normal(self.nbins)
        noisy.sigma_ghat = self.deltag
        noisy.beta_s = beta_s
        noisy.beta_s2 = beta_s**2

        return noisy

####

class HSTBinnerWrapper(HSTBinning):

    def configure(self, config):
        self.parent = config['hstbinning']

    def __call__(self, *args, **kwds):
        return self.parent.doBinning(*args, **kwds)

####

class HSTBinNoiserWrapper(HSTBinning):

    def configure(self, config):
        self.parent = config['hstbinning']

    def __call__(self, *args, **kwds):
        return self.parent.addNoise(*args, **kwds)
