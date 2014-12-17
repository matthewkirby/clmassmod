#####################
# Exploration of noise issues in perfect NFW fits
#####################


import numpy as np
import nfwfit, nfwutils
import nfwmodeltools as tools
import fitmodel
import scipy.misc
import basicMassCon




######################



def createPerfectProfile(m200, c, zcluster, r_mpc, beta_s):

    rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2


    r_scale = nfwutils.rscaleConstM(m200, c, zcluster, 200)

    nfw_shear_inf = tools.NFWShear(r_mpc, c, r_scale, rho_c_over_sigma_c)
    nfw_kappa_inf = tools.NFWKappa(r_mpc, c, r_scale, rho_c_over_sigma_c)
        
    g = beta_s*nfw_shear_inf / (1 - (beta_s*nfw_kappa_inf))

    return g



#########################

def createClusterSet(config, mtrues, zcluster, r_mpc, beta_s, noise, concens = None):

    nclusters = len(mtrues)

    model = nfwfit.buildModel(config)

    if concens is None:

        try:
            mcrelation = model.massconRelation
        except AttributeError:
            mcrelation = nfwfit.buildObject('basicMassCon', 'Duffy', config)


    shearprofiles = []
    shearerr = (noise/np.sqrt(r_mpc))*np.ones_like(r_mpc)


    for i in range(nclusters):
        
        m200 = mtrues[i]

        if concens is None:
            c200 = mcrelation(m200*nfwutils.global_cosmology.h, zcluster, 200)
        else:
            c200 = concens[i]

        gtrue = createPerfectProfile(m200, c200,
                                     zcluster, r_mpc, beta_s)

        shearprofiles.append(gtrue + shearerr*np.random.standard_normal(size=len(r_mpc)))


    return shearprofiles, shearerr


##########################

def fitClusterSet(config, zcluster, r_mpc, beta_s, shearprofiles, shearerr):

    fitter = nfwfit.buildFitter(config)

    fitMasses = np.zeros(len(shearprofiles))
    fitErrs = np.zeros(len(shearprofiles))
    mask = np.ones(len(shearprofiles))

    for i in range(len(shearprofiles)):

        fitres = fitter.minChisqMethod(r_mpc, shearprofiles[i], shearerr, beta_s, beta_s**2, zcluster)
        if fitres is None:
            fitres = fitter.minChisqMethod(r_mpc, shearprofiles[i], shearerr, beta_s, beta_s**2, 
                                           zcluster, useSimplex=True)

        if fitres is None:
            mask[i] = 0
        else:

            

            fitMasses[i] = fitres[0]['m200']*fitter.model.massScale
#            fitMasses[i] = nfwutils.Mdelta(fitres[0]['rscale'], 4.0, zcluster, 200.)
            asymerrs = fitres[1]['m200']
            try:
                fitErrs[i] = fitter.model.massScale*(asymerrs[1] - asymerrs[0])/2.
            except TypeError:
                fitErrs[i] = asymerrs
#
#    mask = (mask == 1)

    return fitMasses, fitErrs, mask

#    return fitMasses, mask


#########################

def bootstrapMean(distro, nboots = 1000):

    means = np.zeros(nboots)
    for i in range(nboots):
        curboot = np.random.randint(0, len(distro), len(distro))
        means[i] = np.mean(distro[curboot])

    return means



    

#########################

def runNoiseScaling(config, nclusters, zcluster, r_mpc, beta_s, noiselevels):

    bias = np.zeros(len(noiselevels))
    biaserr = np.zeros(len(noiselevels))

    mfitsets = []
    merrsets = []
    masksets = []

    mfitsets2 = []
    merrsets2 = []
    masksets2 = []

    mtrues = 10**np.random.uniform(15.8, 16., size=nclusters)

    for i in range(len(noiselevels)):

        shearprofiles, shearerr = createClusterSet(config, mtrues, zcluster, r_mpc, beta_s, noiselevels[i])

        fitMasses, mask  = fitClusterSet(config, zcluster, r_mpc, beta_s, shearprofiles, shearerr)
        mfitsets.append(fitMasses)
#        merrsets.append(fitErrs)
        masksets.append(mask)


#        fitMasses2, fitErrs2, mask2 = fitClusterSet(config, zcluster, r_mpc, beta_s, shearprofiles, shearerr/10.)
#        mfitsets2.append(fitMasses2)
#        merrsets2.append(fitErrs2)
#        masksets2.append(mask2)
#    return mtrues, mfitsets, merrsets, masksets, mfitsets2, merrsets2, masksets2

#    return mtrues, mfitsets, merrsets, masksets
    return mtrues, mfitsets, masksets

    
##########################


def calcNoiseBias(zcluster=0.5, concen=4.):

    m200s = np.arange(1, 50)
    r_mpc = np.arange(0.1, 3.0, 0.1)

    duffy = basicMassCon.Duffy(None)

    background = 1.0
    beta_s = nfwutils.global_cosmology.beta_s([background], zcluster)
    beta_s2 = beta_s**2
    beta_cor = beta_s2/beta_s

    rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2

    def reducedShear(scaledm200):
        m200 = 1e14*scaledm200
        c = duffy(m200, zcluster)
        rs = nfwutils.rscaleConstM(m200, c, zcluster, 200)

        gamma = tools.NFWShear(r_mpc, c, rs, rho_c_over_sigma_c)
        kappa = tools.NFWKappa(r_mpc, c, rs, rho_c_over_sigma_c)

        g = beta_s*gamma/(1-beta_cor*kappa)

        return g

    shearprime = np.zeros((len(m200s), len(r_mpc)))
    shearprimeprime = np.zeros((len(m200s), len(r_mpc)))
    for i, m200 in enumerate(m200s):

        shearprime[i,:] = scipy.misc.derivative(reducedShear, m200, dx = .01, order=5)
        shearprimeprime[i,:] = scipy.misc.derivative(reducedShear, m200, dx = .01, n=2, order=7)
            
    return shearprime, shearprimeprime, -0.5*shearprimeprime/(shearprime**3)


#########################

class Likelihood(object):

    def __init__(self, config):

        self.config = config
        self.model = nfwfit.buildModel(config)

    ###

    
    def bindData(self, r_mpc, g, gerr, beta_s, zcluster):

        self.r_mpc = r_mpc
        self.g = g
        self.gerr = gerr
        self.beta_s = beta_s
        self.beta_s2 = self.beta_s**2
        self.zcluster = zcluster

        self.rho_c = np.float64(nfwutils.global_cosmology.rho_crit(zcluster))
        self.rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2


    ###

    def __call__(self, m200, c200 = None):

        if c200 is None:
            c200 = self.model.massconRelation(m200*nfwutils.global_cosmology.h, self.zcluster, 200)

        return tools.shearprofile_like(m200,
                                       c200,
                                       self.r_mpc,
                                       self.g,
                                       self.gerr,
                                       self.beta_s,
                                       self.beta_s2,
                                       self.rho_c,
                                       self.rho_c_over_sigma_c)


        

#############################


def calcMassUncert(config, r_mpc, beta_s, zcluster, shearprofiles, shearerr,     
                   testmasses = np.arange(-5.02e15, 5e15, 5e13)):

    likelihood = Likelihood(config)

    nclusters = len(shearprofiles)

    testconcens = np.linspace(1.1, 19.9, 100)
    CGrid, MGrid = np.meshgrid(testconcens, testmasses)
    logprobgrid = np.zeros_like(CGrid)
    probgrid = np.zeros_like(logprobgrid)
    margprobgrid = np.zeros_like(testmasses)

    massprobs = np.zeros((nclusters, len(testmasses)))

    masserrs = np.zeros(nclusters)

    for curcluster in range(nclusters):

        likelihood.bindData(r_mpc, shearprofiles[curcluster], shearerr, beta_s, zcluster)

        for i in range(logprobgrid.shape[0]):
            for j in range(logprobgrid.shape[1]):
                logprobgrid[i,j] = likelihood(MGrid[i,j], CGrid[i,j])

        probgrid[:] = np.exp(logprobgrid - np.max(logprobgrid))
        probgrid /= np.sum(probgrid)
        margprobgrid[:] = np.sum(probgrid, axis=1)
        
        massprobs[curcluster,:] = margprobgrid

        avemass = np.sum(margprobgrid*testmasses)
        masserr = np.sqrt(np.sum(margprobgrid*(testmasses - avemass)**2))
        masserrs[curcluster] = masserr


    return masserrs, massprobs

        
        



    
